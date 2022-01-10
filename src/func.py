def get_unique_center_for_cluster(cluster,member=0):
    centers=shapes.loc[cluster[('Alt','Alt1','ID_CENT')]]
    if type(centers)==pd.core.frame.DataFrame:
        if member==0:
            centers=centers.iloc[1]
        if member!=0:
            centers=centers[centers['MEM_MATCH_ID']==member['MEM_MATCH_ID']]
    elif type(centers)==pd.core.series.Series:
        pass
    else:
        raise TypeError
    
    return(centers)

def get_cluster_for_member(clusters,member):
    cluster_id=member[('All','MEM_MATCH_ID')]
    cluster=clusters.loc[cluster_id]
    return(cluster)

def sort2n(x,y):
    """
    Sorts and matches two arrays of object ids where x is unique and y is not (in DES this is coadd_objects_id).
    Slower than sort2().
    """
    
    xsort = np.argsort(x)
    ysort = np.argsort(y)
    i_yx = np.sort(y[np.in1d(y, x, assume_unique=False)])
    i_x = xsort[x[xsort].searchsorted(i_yx)]
    i_y = ysort[y[ysort].searchsorted(i_yx)]
    
    return i_x, i_y


# Drop clusters with no center
def drop_poor_centers(clusters_cat):
    print("The number of clusters before masking", len(clusters_cat))
    center_id=clusters_cat[('Alt', 'Alt1', 'ID_CENT')]
    shape_index=shapes.index
    valid_clusters=np.isin(center_id,shape_index)
    cluster_masked=clusters_cat[valid_clusters]
    print("The number of clusters with no center shape data {}".format(np.sum(~valid_clusters)))
    print("The number of clusters after masking", len(cluster_masked))
    return(cluster_masked)

# clusters=drop_poor_centers(clusters)

# Drop members with no cluster
def drop_poor_members(clusters_cat,shapes_cat):
    print("The number of members before masking: ",len(shapes_cat))
    valid_clusters_id=clusters_cat.index.to_numpy()
    shapes_match_id=shapes_cat[('All','MEM_MATCH_ID')].to_numpy()
    valid_members=np.isin(shapes_match_id,valid_clusters_id)
    print("The number of members without a cluster: ",np.sum(~valid_members))
    members_masked=shapes_cat[valid_members]
    print("Number of galaxies after masking: {}".format(len(members_masked)))
    return(members_masked)
    
# shapes=drop_poor_members(clusters,shapes)

