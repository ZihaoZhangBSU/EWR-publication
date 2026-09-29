"""Execute unchanged notebook code with explicit filesystem-setting overrides.

Only module-level assignments to named settings are replaced. Function bodies,
formulas, thresholds, masks and plot settings are executed exactly as supplied.
"""
from __future__ import annotations
import argparse
import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def execute(path, overrides=None, font=None):
    path = Path(path).resolve()
    overrides = dict(overrides or {})
    overrides.setdefault('NOTEBOOK_DIR',path.parent)
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    if font:
        font_manager.fontManager.addfont(str(Path(font).resolve()))
    sys.path.insert(0, str(ROOT/'workflow/src'))
    sys.path.insert(0, str(path.parent))
    # Notebook cells may display DataFrames; do not require an IPython kernel.
    namespace = {'__name__':'__main__', '__file__':str(path), 'display':print,
                 '_paper_settings':overrides}
    used=set()
    notebook=json.loads(path.read_text(encoding='utf-8'))
    for i, cell in enumerate(notebook['cells']):
        if cell['cell_type'] != 'code': continue
        source=''.join(cell['source'])
        tree=ast.parse(source)
        for node in tree.body:
            if isinstance(node, ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0], ast.Name):
                name=node.targets[0].id
                if name in overrides:
                    node.value=ast.Subscript(value=ast.Name(id='_paper_settings',ctx=ast.Load()),slice=ast.Constant(name),ctx=ast.Load())
                    used.add(name)
        ast.fix_missing_locations(tree)
        print(f'Running {path.name}, cell {i+1}',flush=True)
        exec(compile(tree,f'{path.name}:cell{i+1}','exec'),namespace)
    if set(overrides)-used:
        print('Settings absent from this notebook:', ', '.join(sorted(set(overrides)-used)))
    return namespace

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebook',type=Path)
    parser.add_argument('--workflow-root',type=Path,help='Existing v0.8.2 directory containing src/ and config/.')
    parser.add_argument('--output-root',type=Path,help='Existing output1 directory containing tables/, cache/ and regridded/. Also receives derived outputs.')
    parser.add_argument('--font',type=Path,default=ROOT/'analysis/wy2026/fonts/NimbusSans-Regular.ttf')
    parser.add_argument('--set-path',nargs=2,action='append',default=[],metavar=('SETTING','PATH'))
    parser.add_argument('--set-json',nargs=2,action='append',default=[],metavar=('SETTING','JSON'))
    args=parser.parse_args()
    overrides={}
    if args.workflow_root:
        root=args.workflow_root.expanduser().resolve()
        if not (root/'src/east_river_workflow').is_dir() or not (root/'config/east_river_config.yaml').is_file():
            parser.error('--workflow-root must contain src/east_river_workflow and config/east_river_config.yaml')
        overrides.update(ROOT=root, REPO_ROOT=root, WORKFLOW_ROOT=root, REPO_CANDIDATES=[root], CONFIG_PATH=root/'config/east_river_config.yaml')
    if args.output_root:
        root=args.output_root.expanduser().resolve()
        if not (root/'tables').is_dir():parser.error('--output-root must contain the existing input tables/')
        overrides.update(OUTPUT_DIR=root,WORKFLOW_OUTPUT_DIR=root,PRODUCTION_OUTPUT_DIR=root)
    for name,value in args.set_path:overrides[name]=Path(value).expanduser().resolve()
    for name,value in args.set_json:overrides[name]=json.loads(value)
    execute(args.notebook,overrides,args.font)

if __name__=='__main__':main()
