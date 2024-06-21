from gview import visual
from siman.calc_manage import smart_structure_read

gb = smart_structure_read('test_structure.POSCAR')
ws = visual(gb)

ws.plot(param=dict(
    bonds_color_scale='Viridis',
    projection=True,
    bonds=True,
    bonds_color=True,
    bonds_length=False,
    cell_vectors=True,
))

ws.fig.show()