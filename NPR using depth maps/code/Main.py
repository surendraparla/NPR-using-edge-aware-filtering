import cv2
import numpy
from InputHandling import take_input
from stylize import normal_stylize
from stylizeBW import stylize_bw
from stylizeLines import stylize_lines

typeI,filename,N,points = take_input()

if(typeI == 1):
    normal_stylize(filename,N,points)
elif(typeI == 2):
    stylize_bw(filename,N,points)
elif(typeI == 3):
    stylize_lines(filename,N,points)
