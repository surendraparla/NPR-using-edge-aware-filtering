# NPR-using-edge-aware-filtering
## About
Non-photorealistic rendering (NPR) is an area of computer
graphics that focuses on enabling a wide variety of expressive
styles for digital art. The input to a two dimensional NPR
system is typically an image or video. The output is a typically
an artistic rendering of that input imagery (for example in a
watercolor, painterly or sketched style).In our project we wish
to accomplish this by using edge aware filters. Edge aware
filters blurs the image without effecting the edges.

## Commands to run the project
```
Image = imread('taj.jpg');
imshow(Image);
Out = cartoonize(Image);
figure, imshow(Out);
```
The additional parameters are mentioned in the code and commented.
The code is well documented and self explanatory.

# NPR using depth maps
## Requirements
Python3 installed with numpy and openCV modules.

## How to run
1) Run main.py
2) specify type of stylization you want.
3) Give name of the input file with extension.Make sure that the image is in same directory as the file main.py
4) Specify number of light sources with their coordinates.
5) an output image with name filename_out.jpg will be created in the folder.
### Example:
```
python main.py
Please enter the number of type of stylization you want:
1)normal stylization
2)black and white
3)line stylization
1
Enter filename with extension :Buddha.jpg
Enter number of light sources :1
give co-ordinates of light source-1 :10,10,1.2
```
