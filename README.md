# gview
Basic set up for crystal structure visualization using plotly.


# installation
`git clone https://github.com/geonda/gview.git && cd gview`

`python setup.py install`

# usage
The code takes [siman](https://github.com/dimonaks/siman.git) structure object as an input. The later supports formats from VASP(POSCAR) to CIF and XYZ. 
```
from siman.calc_manage import smart_structure_read
atoms=smart_structure_read('test.cif')
```
Finally to run the code:
```
from gview import visual

visual_instance=visual(atoms)
visual_instance.plot()
visual_instance.fig.show()
```
![img](https://github.com/geonda/gview/blob/main/example/fig1a.png)

# comments 
The plot method of visual object takes dictionary with parameters allowing basic customization. 
```
.plot(param=dict(bonds=True,
                bonds_length=True,
                bonds_color_scale='Viridis',
                scatter=False,
                bonds=True,
                bonds_color=True,
                bonds_length=False,
                cell_vectors=False,
                resolution=50))
```

Once hovered over atoms or bonds it shows related info. One can also choose which bonds to show specifing option `bond_elements=['Ni', 'O']`.

![](https://github.com/geonda/gview/blob/main/example/fig2.png)

One can specify 2d vs 3d view by turning on projection flag. 
![](https://github.com/geonda/gview/blob/main/example/fig3.png)
Inherited from ploty `.fig` attribute can be use for futher customization. 

