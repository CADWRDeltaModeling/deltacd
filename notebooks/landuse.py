#%%
import pandas as pd
import geopandas as gpd
# viz imports
import geoviews as gv
import hvplot.pandas
import holoviews as hv
from holoviews import opts
hv.extension('bokeh')
#
import panel as pn
pn.extension()
#%%
# Read the land use data
landuse = pd.read_csv('../examples/inputs/landuse_dsm2.csv')
landuse.head()
#%%
# Load geojson for land use data
subareas = gpd.read_file('./maps/detaw_168_subareas.geojson').to_crs(epsg=26910)
subareas.head()
subareas.hvplot.polygons(geo=True, crs=26910, tiles='CartoLight', width=500, height=700)
#%%
year=2015
landuse_yr = landuse[landuse.DATE==year]
subarea_yr = subareas.merge(landuse_yr, right_on='area_id', left_on='NEW_SUB')
#%%
crop_type='WS'
subarea_yr[f'{crop_type}_PCT']=subarea_yr[crop_type]/subarea_yr['ACRES']
map = pn.Row(subarea_yr.hvplot.polygons(geo=True, tiles='CartoLight', crs=26910, hover_cols='all', c=crop_type+'_PCT', cmap='viridis',
                                         alpha=0.5, legend=True, rasterize=True, width=500, height=700),)

map

# %%
