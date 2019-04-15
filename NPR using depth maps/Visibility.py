##defining the visibility function
def visibility(p,q,depth):
    x1,y1,z1 = p
    x2,y2,z2 = q
    for x in range(int(min(x1,x2))+1,int(max(x1,x2))):
        #calculate the corresponding y_coordinate.
        slope_y = (y2-y1)/(x2-x1)
        y = y1+slope_y*(x-x1)
        y = int(y)
        height = depth[x][y]

        #Calculate the corresponding z coordinate.
        slope_z = (z2-z1)/(x2-x1)
        z = z1+slope_z*(x-x1)
        #print(height,z)
        if(height > z):
            return False
    return True
