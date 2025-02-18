import math
from collections.abc import Iterable, Sequence
from copy import deepcopy
from pathlib import Path

import openmm
import openmm.app
import openmm.unit
import rdkit
from openff.toolkit import ForceField, Molecule, Quantity, Topology
from openff.toolkit.utils.exceptions import InvalidAtomMetadataError
from openff.toolkit.utils.toolkits import NAGLToolkitWrapper
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from openff.interchange import Interchange
from openff.interchange.components.potentials import Potential
from openff.interchange.exceptions import NonIntegralMoleculeChargeError
from openff.interchange.models import (
    LibraryChargeTopologyKey,
    PotentialKey,
    SingleAtomChargeTopologyKey,
)
from openff.pablo._utils import draw_molecule

__all__ = [
    "draw_molecule",
    "nglview_show_openmm",
    "react",
    "get_openmm_total_charge",
]


def nglview_show_openmm(
    topology, positions: str | Path | Quantity, image_molecules=False
):
    import mdtraj
    import nglview
    import numpy as np
    from openff.units import ensure_quantity

    top = mdtraj.Topology.from_openmm(topology)

    if isinstance(positions, str) or isinstance(positions, Path):
        traj = mdtraj.load(positions, top=top)
        if image_molecules:
            traj.image_molecules(inplace=True)
    else:
        positions = ensure_quantity(positions, "openmm").value_in_unit(
            openmm.unit.nanometer
        )
        xyz = np.asarray([positions])
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            l1, l2, l3, alpha, beta, gamma = (
                mdtraj.utils.box_vectors_to_lengths_and_angles(
                    *np.asarray(box_vectors.value_in_unit(openmm.unit.nanometer))
                )
            )
            unitcell_angles, unitcell_lengths = [alpha, beta, gamma], [l1, l2, l3]
        else:
            unitcell_angles, unitcell_lengths = None, None
        traj = mdtraj.Trajectory(
            xyz, top, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles
        )
    widget = nglview.show_mdtraj(traj)
    widget.clear_representations()
    widget.add_cartoon()
    widget.add_line(opacity=0.5, crossSize=1.0)
    return widget


def get_charge_sum(
    interchange: Interchange,
    topology_indices: Sequence[int],
) -> Quantity:
    charges = {
        key.atom_indices[0]: value.m
        for key, value in interchange["Electrostatics"].charges.items()
    }

    return Quantity(
        sum([charges[index] for index in topology_indices]),
        "elementary_charge",
    )


def smear_charges(
    interchange: Interchange,
    topology_indices: Sequence[int],
) -> Interchange:
    total_formal_charge_of_topology_indices = sum(
        [interchange.topology.atom(i).formal_charge for i in topology_indices]
    )

    initial_charge_sum = get_charge_sum(interchange, topology_indices)

    charge_to_smear = initial_charge_sum - total_formal_charge_of_topology_indices

    per_atom_difference = charge_to_smear / len(topology_indices)

    interchange["Electrostatics"]._charges_cached = False

    for index in topology_indices:
        topology_key = SingleAtomChargeTopologyKey(this_atom_index=index)
        potential_key = interchange["Electrostatics"].key_map[topology_key]

        interchange["Electrostatics"].potentials[potential_key].parameters[
            "charge"
        ] -= per_atom_difference

    new_charge_sum = get_charge_sum(interchange, topology_indices)

    assert math.isclose(
        (initial_charge_sum - new_charge_sum).m,
        charge_to_smear.m,
    )

    assert math.isclose(
        new_charge_sum.m,
        total_formal_charge_of_topology_indices.m,
        abs_tol=1e-10,
        rel_tol=0,
    )

    return interchange


def get_openmm_total_charge(system: openmm.System) -> float:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return sum(
                [
                    force.getParticleParameters(index)[0]._value
                    for index in range(force.getNumParticles())
                ]
            )


def parametrize_with_nagl(
    force_field: ForceField,
    topology: Topology,
    nagl_method: str = "openff-gnn-am1bcc-0.1.0-rc.3.pt",
    allow_nonintegral_charges: bool = False,
) -> Interchange:
    print("adding dummy charges to force field ...")
    ff = deepcopy(force_field)
    # Add a dummy 0.0 library charge at the _top_ so it's only used as a last resort
    ff["LibraryCharges"].add_parameter(
        parameter_kwargs={
            "smirks": "[*:1]",
            "charge1": Quantity(0.0, "elementary_charge"),
            "name": "dummy",
        },
        before=0,  # "[#3+1:1]",
    )

    print("making Interchange ...")
    interchange: Interchange = ff.create_interchange(
        topology,
        allow_nonintegral_charges=True,
    )

    # Remove any assigned charges from Interchange so that assigned charges
    # only come from NAGL
    for molecule in interchange.topology.molecules:
        molecule._partial_charges = None

    potential_keys_to_remove = list()
    topology_keys_to_remove = list()

    nagl_indices = tuple(
        key.atom_indices[0]
        for key, val in interchange["Electrostatics"].key_map.items()
        if val.id == "[*:1]"
    )

    print("replacing dummy charges with NAGL charges ... ")
    for key, charge in interchange["Electrostatics"].charges.items():
        if key.atom_indices[0] in nagl_indices:
            # only modify charges where dummy placeholder of 0.0 was assigned from "[*:1]" parameter
            index = key.atom_indices[0]

            # If we haven't seen this molecule before, we need to calculate its
            # partial charges
            atom = interchange.topology.atom(index)
            molecule = atom.molecule
            index_in_molecule = atom.molecule_atom_index
            if molecule.partial_charges is None:
                print(
                    f"assigning graph charges to {molecule.name or molecule.hill_formula} ..."
                )
                molecule.assign_partial_charges(
                    partial_charge_method=nagl_method,
                    toolkit_registry=NAGLToolkitWrapper(),
                )
                print("continuing with dummy charge replacement ...")

            # must make new "single atom" topology key for each atom since the current
            # 1:many representation from the dummy charge is no longer valid
            new_potential_key = PotentialKey(
                id="inserted_graph_charges",
                associated_handler="molecules_with_preset_charges",
                mult=index,
            )
            # Place this atom's NAGL charge in its own potential
            new_potential = Potential(
                parameters={"charge": molecule.partial_charges[index_in_molecule]}
            )

            # Add the new potential and key to the interchange
            interchange["Electrostatics"].key_map[
                SingleAtomChargeTopologyKey(this_atom_index=index)
            ] = new_potential_key
            interchange["Electrostatics"].potentials.update(
                {new_potential_key: new_potential}
            )

            # remove the keys associated with the dummy library charge
            potential_keys_to_remove.append(
                interchange["Electrostatics"].key_map[
                    LibraryChargeTopologyKey(this_atom_index=index)
                ]
            )

            topology_keys_to_remove.append(
                LibraryChargeTopologyKey(this_atom_index=index)
            )

    for key_to_remove in topology_keys_to_remove:
        interchange["Electrostatics"].key_map.pop(key_to_remove)

    interchange["Electrostatics"]._charges_cached = False

    interchange = smear_charges(
        interchange=interchange,
        topology_indices=nagl_indices,
    )

    if not allow_nonintegral_charges:
        total_formal_charge = sum(
            atom.formal_charge for atom in interchange.topology.atoms
        ).m_as("elementary_charge")

        net_charge = get_charge_sum(
            interchange=interchange,
            topology_indices=range(interchange.topology.n_atoms),
        ).m_as("elementary_charge")

        if abs(total_formal_charge - net_charge) > 0.01:
            raise NonIntegralMoleculeChargeError(
                f"Interchange has a net charge of {net_charge} compared to a"
                + f" total formal charge  of {total_formal_charge}.",
            )

    return interchange


def react(
    reactants: Sequence[Molecule],
    reaction_smarts: str,
) -> Iterable[tuple[Molecule, ...]]:
    # Convert reactants to rdmol, storing metadata as properties
    # Need to preserve metadata so we can identify leaving atoms and synonyms
    reactant_rdmols = [reactant.to_rdkit() for reactant in reactants]
    for reactant_rdmol, reactant_offmol in zip(reactant_rdmols, reactants):
        for reactant_rdatom, reactant_offatom in zip(
            reactant_rdmol.GetAtoms(), reactant_offmol.atoms
        ):
            for key, value in reactant_offatom.metadata.items():
                if isinstance(value, bool):
                    reactant_rdatom.SetBoolProp(key, value)
                elif isinstance(value, int):
                    reactant_rdatom.SetIntProp(key, value)
                elif isinstance(value, float):
                    reactant_rdatom.SetDoubleProp(key, value)
                else:
                    reactant_rdatom.SetProp(key, str(value))

    # Prepare the reaction
    rxn = ReactionFromSmarts(reaction_smarts)
    product_rdmols = rxn.RunReactants(reactant_rdmols)

    # Get map from reaction SMARTS atom mappings to the equivalent OFF atom
    map_to_offatom = {}
    for reactant_rdmol, reactant_offmol in zip(reactant_rdmols, reactants):
        assert rxn.IsMoleculeReactant(reactant_rdmol)
        for reactant_template in rxn.GetReactants():
            map_to_offatom.update(
                {
                    reactant_template.GetAtomWithIdx(query).GetProp(
                        "molAtomMapNumber"
                    ): reactant_offmol.atom(match)
                    for query, match in enumerate(
                        reactant_rdmol.GetSubstructMatch(reactant_template)
                    )
                    if reactant_template.GetAtomWithIdx(query).HasProp(
                        "molAtomMapNumber"
                    )
                }
            )

    # Process and yield the products
    for products in product_rdmols:
        # Skip products that cannot be sanitized
        try:
            for product in products:
                product.UpdatePropertyCache()
                rdkit.Chem.SanitizeMol(product)
        except rdkit.Chem.rdchem.MolSanitizeException:
            continue

        product_offmols = [Molecule.from_rdkit(product) for product in products]

        # Fix metadata of products
        for product_template in rxn.GetProducts():
            for product_rdmol, product_offmol in zip(products, product_offmols):
                # Go over atoms changed in the reaction and fix their metadata (rdkit often loses it)
                for product_idx, product_template_idx in enumerate(
                    product_rdmol.GetSubstructMatch(product_template)
                ):
                    product_rdatom = product_rdmol.GetAtomWithIdx(product_idx)
                    product_offatom = product_offmol.atom(product_idx)

                    if product_rdatom.HasProp("old_mapno"):
                        rxn_map = product_rdatom.GetProp("old_mapno")
                        reactant_offatom = map_to_offatom[rxn_map]
                        product_offatom.metadata.update(reactant_offatom.metadata)
                        if "leaving_atom" in product_offatom.metadata:
                            product_offatom.metadata["leaving_atom"] = False
                        if "substructure_atom" in product_offatom.metadata:
                            product_offatom.metadata["substructure_atom"] = True
                        product_offatom.name = reactant_offatom.name

                # Copy the props back to the metadata
                for product_rdatom, product_offatom in zip(
                    product_rdmol.GetAtoms(), product_offmol.atoms
                ):
                    for key, value in product_rdatom.GetPropsAsDict().items():
                        try:
                            product_offatom.metadata[key] = value
                        except InvalidAtomMetadataError:
                            pass

        yield tuple(product_offmols)
