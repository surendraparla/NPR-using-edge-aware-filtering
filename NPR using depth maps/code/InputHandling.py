import numpy
def take_input():
    #Taking input from the user.
    message = 'Please enter the number of type of stylization you want:\n'+'1)normal stylization\n'+'2)black and white\n'+'3)line stylization \n'
    typeI = int(input(message))
    filename = str(input('Enter filename with extension :'))
    N = int(input('Enter number of light sources :'))
    points = numpy.zeros((N,3))
    for i in range(N):
        message = 'give co-ordinates of light source-'+str(i+1)+' :'
        points[i] = [float(x) for x in input(message).split(',')]

    return typeI,filename,N,points
