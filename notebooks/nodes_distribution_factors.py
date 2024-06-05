
# %%
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
# %%
nodes_area_file = "./maps/nodes_in_areas.geojson"
old_div_file = "../examples/inputs/diversion_factors_dsm2_mss.csv"
new_div_file = "../examples/inputs/diversion_factors_dsm2_mss_equalized_factor.csv"

# %%
nodes_areas_all_df = gpd.read_file(nodes_area_file).to_crs(epsg=26910)
old_div_df = pd.read_csv(old_div_file, index_col=False)
# %%
# save split area (BBID) in a different datafram
new_split_area_df = old_div_df.where(old_div_df['area_id']>168)
# select the 168 subareas
old_div_168_df = old_div_df.where(old_div_df['area_id']<=168)
# drop BBID
old_div_168_df.dropna(inplace=True)
new_split_area_df.dropna(inplace=True)


# %%
nodes_areas_df = nodes_areas_all_df[['id','NEW_SUB']]
nodes_areas_df.rename(columns={'id': 'node','NEW_SUB': 'area_id'}, inplace=True)
nodes_areas_df.sort_values(by=['area_id','node'],inplace=True)

# The original csv doesn't read node as int becase it has BBID in it. Hence, need to force change type to it.
old_div_168_df['node']=old_div_168_df['node'].astype('Int64')
# %%

new_dist_df = pd.merge(nodes_areas_df,old_div_168_df,on=['area_id','node'],how='outer')

# %%
# calculate the factor
new_dist_value_counts = new_dist_df['area_id'].value_counts()
new_dist_value_counts.sort_index(inplace=True)
new_factor = 1/new_dist_value_counts
# repeat the factor to match the length of the dataframe
new_factor_df = pd.DataFrame(new_factor.repeat(new_dist_value_counts.values).values,columns=['factor'])
# replace column with new factor
new_dist_df['factor'] = new_factor_df['factor']

# reorder the columns
new_dist_df = new_dist_df[['area_id','node','factor']]
# concat the split area dataframe
new_dist_df = pd.concat([new_dist_df,new_split_area_df])
# %%
# write to csv
new_dist_df.to_csv(new_div_file,float_format="%.4f",index=False)
