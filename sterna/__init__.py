bl_info = {
    'name': 'Sterna',
    'author': 'Antti Elonen',
    'category': 'All',
    'version': (1, 0, 20200203),
    'blender': (2, 82, 0)
}

try:
    import importlib
    importlib.reload(blenderEdgeGroups)
    importlib.reload(sternaMain)
    importlib.reload(sternaEdit)
    importlib.reload(sternaExport)
    importlib.reload(sternaImport)
    importlib.reload(propDefs)
    importlib.reload(sternaStrand)
except:
    from . import blenderEdgeGroups, sternaMain, sternaEdit, sternaExport, sternaImport, propDefs, sternaStrand

import sys, os, bpy

def register():
    blenderEdgeGroups.register()
    sternaMain.register()
    sternaEdit.register()
    sternaExport.register()
    sternaImport.register()

def unregister():
    blenderEdgeGroups.unregister()
    sternaMain.unregister()
    sternaEdit.unregister()
    sternaExport.unregister()
    sternaImport.unregister()

if __name__ == "__main__":
    register()
