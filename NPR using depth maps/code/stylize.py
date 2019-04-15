import cv2
import numpy
from Visibility import visibility

def normal_stylize(filename,N,points):
    # Creating the L-channel
    img = cv2.imread(filename)
    cv2.imshow(filename[:-4]+'_input',img)
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

    #Creating new image based on visibility.
    d_img = numpy.zeros((N,img_L.shape[0],img_L.shape[1]),dtype = img_L.dtype)
    for n in range(N):
        for i in range(img_L.shape[0]):
            for j in range(img_L.shape[1]):
                q = (i,j,depth[i][j])
                if(visibility(points[n],q,depth)):
                    d_img[n][i][j] = img_L[i][j]

    f_d_img = numpy.zeros(img_L.shape,dtype = 'float64')
    for i in range(N):
        f_d_img += d_img[i]/N 

    f_img = numpy.zeros(img_Lab.shape,dtype = img_Lab.dtype)
    f_img[:,:,0] = f_d_img.astype(img_L.dtype)
    f_img[:,:,1] = img_A[:,:]
    f_img[:,:,2] = img_B[:,:]
    f_img = cv2.cvtColor(f_img,cv2.COLOR_Lab2BGR)
    cv2.imshow(filename[:-4]+'_final',f_img)
    cv2.imwrite(filename[:-4]+'_out.jpg',f_img)





