#%%
import xarray as xr
import pandas as pd
import geopandas as gpd
from functools import partial
# viz imports
import geoviews as gv
import hvplot.pandas
import hvplot.xarray
import holoviews as hv
from holoviews import opts, dim
hv.extension('bokeh')
#
import panel as pn
pn.extension()

# %%
# Create a map of the nodes with sizes based on the diversion amount 
ncfile = '../examples/output/dcd_dsm2_mss1.nc'
nodes_file = './maps/dsm2_nodes_8_2.geojson'
# %%
# read the nc data
xrdata = xr.open_dataset(ncfile)
# %%
dflist = {}
for var in ['diversion', 'seepage', 'drainage']:
    monthly = xrdata[var].resample(time='M').mean()
    monthly = monthly.to_dataframe().reset_index()
    df=monthly.pivot(index='time', columns='node', values=var)
    df = df.reset_index()
    df.melt(id_vars='time', value_vars=df.columns)
    # BBID mapped on to node 94 (There was no value from the diversion data for that node so substituting with BBID)
    df = df.rename(columns={'BBID':'94'})
    dflist[var] = df
# %%
gdf = gpd.read_file(nodes_file).to_crs(epsg=26910)
# %%
import numpy as np
def build_map(time, df=None, var=''):
    if var == 'diversion':
        color='blue'
    elif var == 'seepage':
        color='green'
    else:
        color='red'
    dft = df[df['time']==time].melt(id_vars='time', value_vars=df.columns)
    dft['node'] = dft['node'].astype(int)
    gdfm = gdf.merge(dft, left_on='id', right_on='node')
    return gdfm.hvplot.points(geo=True, crs=26910, hover_cols='all', alpha=0.35, width=300, tiles='CartoLight', c=color,
                              responsive=False, legend=True).opts(opts.Points(size=0.5*dim('value'))).opts(title=var.upper()+' '+pd.to_datetime(time).strftime('%Y-%m'))
#
# %%
times = dflist['diversion']['time'].unique()
#%%
dmaplist = {var:hv.DynamicMap(partial(build_map, df=dflist[var], var=var), kdims='time').redim.values(time=times) for var in dflist.keys()}
maplist = pn.Row(dmaplist['diversion'], dmaplist['drainage'], dmaplist['seepage']).servable()
maplist.show()
# %%
