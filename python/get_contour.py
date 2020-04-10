
def get_contour_data(cont_obj, i):
    # path object
    coll = cont_obj.collections[i]
    try:
        p = coll.get_paths()[0]
    except:
        print 'no path found !'
        return [],[]
    v = p.vertices
    x = v[:,0]
    y = v[:,1]
    return x,y

def write_contour(name, cont_obj, i, var_names = ['x','y']):
    # write contour level i as a data file woth filename name
    level_val = cont_obj.levels[i]
    o = open(name, 'w')
    # write header
    o.write( "#\\ level = {0}\n".format(level_val) )
    o.write ("#! {0}[f,0]/ {1}[f,1]/ \n".format( var_names[0], var_names[1] ) )
    # get the data
    x,y = get_contour_data(cont_obj, i)
    for i,xx in enumerate(x):
        yy = y[i]
        o.write("{0} {1}\n".format(xx, yy) )
    o.close()

#done
    



