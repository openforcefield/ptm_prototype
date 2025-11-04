from pathlib import Path
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

reactants_smarts = [
    "[*]-[S]-[H]",
    "[*]-[N]1-[C](=[O])-[CH]=[CH]-[C](=[O])-1",
]
products_smarts = [
    "[*]-[N]1-[C](=[O])-[CH2]-[CH](-[S]-[*])-[C](=[O])-1",
]
thiol_maleimide_click_smarts = (
    ".".join(reactants_smarts) + ">>" + ".".join(products_smarts)
)

rxn = ReactionFromSmarts(thiol_maleimide_click_smarts)
d2d = Draw.MolDraw2DCairo(-1, 150)
dopts = d2d.drawOptions()
dopts.dummiesAreAttachments = True
dopts.padding = 0.0
d2d.DrawReaction(rxn)
Path("reaction.png").write_bytes(d2d.GetDrawingText())
