# This script finds the nodes that are within a certain distance of the border of the areas
# It uses geopandas and shapely libraries to perform spatial operations.
# The script reads nodes and areas data from GeoJSON  files, creates a buffer around the nodes,
# and then finds the buffered nodes that intersect with the areas.
# The resulting nodes are saved as a new GeoJSON file.

# %%
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon

# %%
nodes_file = "./maps/dsm2_nodes_8_2.geojson"
areas_file = "./maps/detaw_168_subareas.geojson"
# %%
nodes = gpd.read_file(nodes_file).to_crs(epsg=26910)
# create a dataframe with geometry as a buffer around the nodes
buffer_size = 1000  # 1000 ft buffer, units same as the crs
buffered_nodes = nodes.copy()
buffered_nodes.geometry = nodes.buffer(buffer_size)

# %%
areas = gpd.read_file(areas_file).to_crs(epsg=26910)
# %%
# find the buffered nodes that intersect with the areas
nodes_in_areas = gpd.sjoin(buffered_nodes, areas, op="intersects")
# %%
# save as geojson
nodes_in_areas.to_file("./maps/nodes_in_areas.geojson", driver="GeoJSON")
# %%
