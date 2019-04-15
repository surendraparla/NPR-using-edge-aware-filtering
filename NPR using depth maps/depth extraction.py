import cv2
import numpy

##defining the visibility function
def visibility(p,q,depth):
    x1,y1,z1 = p
    x2,y2,z2 = q
    for x in range(min(x1,x2)+1,max(x1,x2)):
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

# Creating the L-channel
img = cv2.imread('buddha.jpg')
cv2.imshow('buddha_input',img)
img_Lab = cv2.cvtColor(img,cv2.COLOR_BGR2Lab)
img_L = img_Lab[:,:,0]
img_A = img_Lab[:,:,1]
img_B = img_Lab[:,:,2]

#Creating Base and depth layers
base = cv2.bilateralFilter(img_L,9,75,75)
detail = cv2.subtract(img_L,base)

#Creating the depth map.
F_b = 0.8
F_d = 0.6

base_N = base/255
detail_N = detail/255
depth = F_b*base_N + F_d*detail_N

#Creating a point source of light
x_co = 10
y_co = 10
h1 = 1.4
p = (x_co,y_co, h1)

#Creating new image based on visibility.

d_img = numpy.zeros(img_L.shape,dtype = img_L.dtype)

for i in range(img_L.shape[0]):
    for j in range(img_L.shape[1]):
        q = (i,j,depth[i][j])
        if(visibility(p,q,depth)):
            d_img[i][j] = img_L[i][j]

f_img = numpy.zeros(img_Lab.shape,dtype = img_Lab.dtype)
f_img[:,:,0] = d_img[:,:]
f_img[:,:,1] = img_A[:,:]
f_img[:,:,2] = img_B[:,:]
f_img = cv2.cvtColor(f_img,cv2.COLOR_Lab2BGR)
print(d_img.dtype)
cv2.imshow('buddha_final',f_img)
cv2.imwrite('buddha_out.jpg',f_img)






