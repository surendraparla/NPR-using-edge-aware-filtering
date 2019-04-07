# NPR-using-edge-aware-filtering
Non-photorealistic rendering (NPR) is an area of computer
graphics that focuses on enabling a wide variety of expressive
styles for digital art. The input to a two dimensional NPR
system is typically an image or video. The output is a typically
an artistic rendering of that input imagery (for example in a
watercolor, painterly or sketched style).In our project we wish
to accomplish this by using edge aware filters. Edge aware
filters blurs the image without effecting the edges.

Commands to run the project

Image = imread('taj.jpg');
imshow(Image);
Out = cartoonize(Image);
figure, imshow(Out);

The additional parameters are mentioned in the code and commented.
The code is well documented and self explanatory.
