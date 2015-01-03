"""
Module Description:
    Graphs single graphs in 3D and initiates 3D interface that allows user to
    rotate, translate, dilate, stretch, and compress the graph.
    
Global variables:
	obj_list
		store standard coordinates of spaces
	rot_list
		store rotated coordinates of spaces
	trans_list
		stores rotating actions
	trans_list_one
		condenses trans_list in the following format:
		[['r',[10,0,10]],['z',0.75],['m',[23,5,-12]] ... ]
	oselected
		int that keeps track of selected object
	sodict
		relates name of object to position in rot_list and onj_list
	axis_dic
		grid dictionary
		{'vertices':[[]],'xp':[],'xn':[],'yp':[],'yn':[],'zp':[],'zn':[]}
		 
Class space:

	rotate_x(self,degrees)
	rotate_y(self,degrees)
	rotate_z(self,degrees)
	rotate(self,x,y,z,X=0,Y=0,Z=0)
	translate(self,x,y,z)
	dilate(self,factor,X=0,Y=0,Z=0)
	compress(self,x,y,z,X=0,Y=0,Z=0)
	reflect(self,direction,X=0,Y=0,Z=0)
		input 'xy', 'yz', or 'xz' for direction parameter
	
	graph(self,screen,Graphturtle)

Module functions:
	axis_dic_ini()
		initiates global axis_dic
	GraphList(Graphstr, Max, Min, increment) --> space
	PenColor(Point) --> str
	ConvertHex(Rcolor, Gcolor, Bcolor) --> str
	TurtleFit(point) --> (x,y)
	trans_list_update(trans_list,trans_list_one)
	axis_display(vlist) --> list of str
	atan2fixed(y,x) --> float
	convex_hull(plist) --> list of (x,y)
	in_poly(xytup, poly)
	org_to_rot(space) 
	main()
		main program
"""

from copy import *
import tkinter
from tkinter import *
import turtle
from math import *
from equation_solver import *
import time

#initiate global variables
DectoHexDic = {1:'1',2:'2',3:'3',4:'4',5:'5',6:'6',7:'7',8:'8',9:'9',10:'A'\
               ,11:'B',12:'C',13:'D',14:'E',15:'F',0:'0'}


class space(list):
    """class used to store a matrix of rectangular coordinates for rotation.
   must be in the form list of list of list of three nums.
   all spaces in the form of:
   [[[x,y,z,color], ... ], ... ]
    """

    def __init__(self):
        """"""

    def rotate_x(self,degrees):
        """rotate object about x-axis for given number of degrees.
        """
        radians=math.radians(degrees)
        for row in self:
            for point in row:
                if point == []:
                    continue
                y = point[1]
                z = point[2]
                Y = y*math.cos(radians)-z*math.sin(radians)
                Z = y*math.sin(radians)+z*math.cos(radians)
                point[1] = Y
                point[2] = Z

    def rotate_y(self,degrees):
        """rotate object about y-axis for given number of degrees.
        """
        radians=math.radians(degrees)
        for row in self:
            for point in row:
                if point == []:
                    continue 
                x = point[0]
                z = point[2]
                X = x*math.cos(radians)+z*math.sin(radians)
                Z = -1*x*math.sin(radians)+z*math.cos(radians)
                point[0] = X
                point[2] = Z
                
    def rotate_z(self,degrees):
        """rotate object about z-axis for given number of degrees.
        """
        radians=math.radians(degrees)
        for row in self:
            for point in row:
                if point == []:
                    continue 
                x = point[0]
                y = point[1]
                X = x*math.cos(radians)-y*math.sin(radians)
                Y = x*math.sin(radians)+y*math.cos(radians)
                point[0] = X
                point[1] = Y
    def rotate(self,x,y,z,X=0,Y=0,Z=0):
        """rotate object about (X,Y,Z) for specified degrees in x,y,and z directions.
        """
        for row in self:
            for point in row:
                if point == []:
                    continue
                point[0] -= X
                point[1] -= Y
                point[2] -= Z
        self.rotate_x(x)
        self.rotate_y(y)
        self.rotate_z(z)
        for row in self:
            for point in row:
                if point == []:
                    continue
                point[0] += X
                point[1] += Y
                point[2] += Z
    def translate(self,x,y,z):
        """translates self x-units in the positive x direction, y-units in the
        positive y direction, and z-units in the positive z direction.
        """
        for row in self:
            for point in row:
                if point == []:
                    continue 
                point[0] += x
                point[1] += y
                point[2] += z

    def dilate(self,factor,X=0,Y=0,Z=0):
        """Dilate the object by factor with respect to origin (X,Y,Z)
        """
        for row in self:
            for point in row:
                if point == []:
                    continue 
                xone = point[0] - X
                yone = point[1] - Y
                zone = point[2] - Z
                xone = xone*factor
                yone = yone*factor
                zone = zone*factor
                point[0] = X + xone
                point[1] = Y + yone
                point[2] = Z + zone
            
    def compress(self,x,y,z,X=0,Y=0,Z=0):
        """compresses self by the specified factors x,y,z about the axes X,Y,Z.
        """
        for row in self:
            for point in row:
                if point == []:
                    continue 
                xone = point[0] - X
                yone = point[1] - Y
                zone = point[2] - Z
                xone = xone*x
                yone = yone*y
                zone = zone*z
                point[0] = X + xone
                point[1] = Y + yone
                point[2] = Z + zone       

    def reflect(self,direction,X=0,Y=0,Z=0):
        """reflect the object over either the xy-plane, xz-plane, or yz-plane
        at specified origin (X,Y,Z).
        Prerequisite: input 'xy', 'yz', or 'xz' for direction parameter.
        """
        if direction == 'xy':
            for row in self:
                for point in row:
                    if point == []:
                        continue 
                    zone = point[2] - Z
                    zone = -zone
                    point[2] = Z + zone
        elif direction == 'xz':
            for row in self:
                for point in row:
                    if point == []:
                        continue 
                    yone = point[1] - Y
                    yone = -yone
                    point[1] = Y + yone
        elif direction == 'yz':
            for row in self:
                for point in row:
                    if point == []:
                        continue 
                    xone = point[0] - X
                    xone = -xone
                    point[0] = X + xone
        
    def graph(self,screen,Graphturtle):
        """Graphs function on turtle."""
        
        Graphturtle.penup()

        Ylinelist = []
        for rowlist in self:
           
            
            YlinelistRow = []
            
            beginning = True
            
            for point in rowlist:
                if point == [] or point[2]>300 or point[2]<-300 or point[0]>300 or point[0]<-300:
                    Graphturtle.penup()
                    
                    YlinelistRow.append([])                                 
                    
                    continue
              
                Graphtup = TurtleFit(point)

                color = point[3]
                
                if beginning:
                    Graphturtle.goto(Graphtup)
                    Graphturtle.pendown()
                    beginning = False
                else:
                    Graphturtle.pendown()
                    Graphturtle.goto(Graphtup)
                
                Graphturtle.color(color)
                YlinelistRow.append([Graphtup,color])
                #Graphturtle.dot()
                Graphturtle.penup()
                
            Graphturtle.penup()
            Ylinelist.append(YlinelistRow)
            screen.update()
        
        #Draw y lines
        length = len(Ylinelist)
        
        
        for i in range(length):
            beginning = True
            for j in range(length):
                
                point = Ylinelist[j][i]
                
                if point != []:
                    if beginning:
                        Graphturtle.pencolor(point[1])
                        Graphturtle.goto(point[0])
                        Graphturtle.pendown()
                        beginning = False
                        continue
                    Graphturtle.pencolor(point[1])
                    Graphturtle.goto(point[0])
            Graphturtle.penup()

            
#initiate list to store standard coordinates of spaces.
obj_list = []
#initiate list to store rotated coordinates of spaces.
rot_list = []
#initiate list to store rotating actions.
trans_list = []
#initiate a veresion of trans_list that condenses information and is
#used in org_to_rot function.
#Translate_list_one uses the following format:
#[['r',[10,0,10]],['z',0.75],['m',[23,5,-12]] ... ]
trans_list_one = []
#oselected is the num that keeps track of which object is selected.
#starts at 0, equivalent to indexing in obj_list and rot_list.
oselected = 0
#sodict relates the name of an object to its position in rot_list and obj_list.
sodict = {}
#axis_dic initiates the grid in format of {'vertices':[[ ... ]],'xp':[],
#'xn':[],'yp':[],'yn':[],'zp':[],'zn':[]}
def axis_dic_ini():
    """initiates grid"""
    ##xgrid red
    ##ygrid green
    ##zgrid blue
    global axis_dic
    axis_dic = {}
    axis_dic['vertices'] = space()
    axis_dic['vertices'].append([[-300,300,300,'#FFFFFF',0]\
                                 ,[-300,-300,300,'#FFFFFF',1]\
                                 ,[-300,-300,-300,'#FFFFFF',2]\
                                 ,[-300,300,-300,'#FFFFFF',3]\
                                 ,[300,300,300,'#FFFFFF',4]\
                                 ,[300,-300,300,'#FFFFFF',5]\
                                 ,[300,-300,-300,'#FFFFFF',6]\
                                 ,[300,300,-300,'#FFFFFF',7]])
    axis_dic['xp'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([300,i,j,'#FFCCCC'])
        axis_dic['xp'].append(row_list)
    axis_dic['xn'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([-300,i,j,'#FFCCCC'])
        axis_dic['xn'].append(row_list)
    axis_dic['yp'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([i,300,j,'#DFFFA5'])
        axis_dic['yp'].append(row_list)
    axis_dic['yn'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([i,-300,j,'#DFFFA5'])
        axis_dic['yn'].append(row_list) 
    axis_dic['zp'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([i,j,300,'#ADD8E6'])
        axis_dic['zp'].append(row_list)
    axis_dic['zn'] = space()
    for i in range(-300,301,40):
        row_list = []
        for j in range(-300,301,40):
            row_list.append([i,j,-300,'#ADD8E6'])
        axis_dic['zn'].append(row_list)

axis_dic_ini()

def GraphList(Graphstr, Max, Min, increment):
    """(str, dic,num) --> space

    Uses Graphstr and given prarmeters to generate a list of 3D coordinates
    for a graph.
    """
    
    Xlist = []
    x = Max
    while x >= Min:
        Xlist.append(x)
        x -= increment
    
    Ylist = []
    y = Max
    while y >= Min:
        Ylist.append(y)
        y -= increment
        
    XYlist = []
    for Y in Ylist:
        rowlist = []
        for X in Xlist:
            rowlist.append((X,Y))
        XYlist.append(rowlist)

    
    returnsurface = space()
    for rowlist in XYlist:
        rowsurfacelist = []
        for XYtup in rowlist:
            try:
                Xvalue = XYtup[0]
                Yvalue = XYtup[1]
                Zvalue = EquationSolver(Graphstr, Xvalue, Yvalue)
                if Zvalue > Max or Zvalue < Min:
                    rowsurfacelist.append([])
                else:
                    color = PenColor([Xvalue,Yvalue,Zvalue])
                    rowsurfacelist.append([Xvalue,Yvalue,Zvalue,color])
            except (ZeroDivisionError, ValueError):
                #When an error is returned, empty str is appended as place holder.
                rowsurfacelist.append([])
        returnsurface.append(rowsurfacelist)
    
    return returnsurface

def PenColor(Point):
    """(list of three num) --> str 

    Return a string in Hex Code format of the color desired for the given
    3D rectangular coordinates.
    The color of a point is determined based on the Z-coordinate
    with a red gradient.

    Prerequisite: The point must be between max and min.
    """   
    Zvalue = Point[2]
    Zdif = 600
    Zpercent = (Zvalue+300)/Zdif
    
    if Zpercent >= 0.5:
        Zpercent = (Zpercent-0.5)*2
        GBcolor = round(220*Zpercent)
        return ConvertHex(255,GBcolor, GBcolor)
    
    elif Zpercent < 0.5:
        Zpercent = (0.5-Zpercent)*2
        Rdecrease = 170*Zpercent
        Rcolor = round(255-Rdecrease)
        return ConvertHex(Rcolor, 0, 0)

def ConvertHex(Rcolor, Gcolor, Bcolor):
    """(int,int,int) --> str

    Convert given RGB tuple (from 0-255) to the a str representing
    Hex Code for that color.

    >>> ConvertHex(255, 0 ,0)
    '#FF0000'

    >>> ConvertHex(100, 100, 150)
    '#646496'
    """
    returnstr = '#'
    returnstr += DectoHexDic[Rcolor//16]
    returnstr += DectoHexDic[Rcolor%16]
    returnstr += DectoHexDic[Gcolor//16]
    returnstr += DectoHexDic[Gcolor%16]
    returnstr += DectoHexDic[Bcolor//16]
    returnstr += DectoHexDic[Bcolor%16]

    return returnstr

def TurtleFit(point):
    """(tuple) -> tuple (x,y)

    Turtle fit function simply takes pre-fitted coordinates and graphs puts
    them on the canvas.
    """
    X = 300+point[0]
    Y = 300+point[2]

    return (X,Y)

def trans_list_update(trans_list,trans_list_one):
    """(list,list) --> Nonetype

    Update trans_list_one to include new information from trans_list.

    Prerequisite: Everytime trans_list_one is updated, one and only one item
    should have been added to the trans_list
    """
    if len(trans_list) == 0:
        return None
    new_trans = trans_list[-1]
    new_trans_type = new_trans[0]
   
    if len(trans_list_one) == 0 or trans_list_one[-1][0] != new_trans_type:
        if new_trans_type == 'r':
            if new_trans == 'rup':
                trans_list_one.append([new_trans_type,[-20,0,0]])
            elif new_trans == 'rdown':
                trans_list_one.append([new_trans_type,[20,0,0]])
            elif new_trans == 'rleft':
                trans_list_one.append([new_trans_type,[0,0,-20]])
            elif new_trans == 'rright':
                trans_list_one.append([new_trans_type,[0,0,20]])
            elif new_trans == 'rtup':
                trans_list_one.append([new_trans_type,[0,-20,0]])
            elif new_trans == 'rtdown':
                trans_list_one.append([new_trans_type,[0,20,0]])
        elif new_trans_type == 'm':
            if new_trans == 'mup':
                trans_list_one.append([new_trans_type,[0,0,20]])
            elif new_trans == 'mdown':
                trans_list_one.append([new_trans_type,[0,0,-20]])
            elif new_trans == 'mleft':
                trans_list_one.append([new_trans_type,[-20,0,0]])
            elif new_trans == 'mright':
                trans_list_one.append([new_trans_type,[20,0,0]])
            elif new_trans == 'min':
                trans_list_one.append([new_trans_type,[0,20,0]])
            elif new_trans == 'mdown':
                trans_list_one.append([new_trans_type,[0,-20,0]])
        elif new_trans_type == 'z':
            if new_trans == 'zout':
                trans_list_one.append([new_trans_type,0.75])
            elif new_trans == 'zin':
                trans_list_one.append([new_trans_type,1.25])
            
    elif trans_list_one[-1][0] == new_trans_type:
        if new_trans_type == 'r':
            if new_trans == 'rup':
                trans_list_one[-1][1][0] -= 20
            elif new_trans == 'rdown':
                trans_list_one[-1][1][0] += 20
            elif new_trans == 'rleft':
                trans_list_one[-1][1][2] -= 20
            elif new_trans == 'rright':
                trans_list_one[-1][1][2] += 20
            elif new_trans == 'rtup':
                trans_list_one[-1][1][1] -= 20
            elif new_trans == 'rtdown':
                trans_list_one[-1][1][1] += 20
        elif new_trans_type == 'm':
            if new_trans == 'mup':
                trans_list_one[-1][1][2] += 20
            elif new_trans == 'mdown':
                trans_list_one[-1][1][2] -= 20
            elif new_trans == 'mleft':
                trans_list_one[-1][1][0] -= 20
            elif new_trans == 'mright':
                trans_list_one[-1][1][0] += 20
            elif new_trans == 'min':
                trans_list_one[-1][1][1] += 20
            elif new_trans == 'mdown':
                trans_list_one[-1][1][1] -= 20
        elif new_trans_type == 'z':
            if new_trans == 'zout':
                trans_list_one[-1][1] = trans_list_one[-1][1]*0.75
            elif new_trans == 'zin':
                trans_list_one[-1][1] = trans_list_one[-1][1]*1.25

def axis_display(vlist):
    """(space) -> list of str

    Return a list of str ('xn','xp','yp','yn','zn','zp') that indicate which
    grid face will be displayed.
    """
    vlist = vlist[0]
    #sort the coordinates by y coordinate (v[0][1]) from big to small
    sortedvlist = [vlist[0]]
    for i in range(1,8):
        j = 0
        length = len(sortedvlist)
        while j < length:
            if vlist[i][1] >= sortedvlist[j][1]:
                sortedvlist.insert(j,vlist[i])
                break
            j += 1
        if j==length:
            sortedvlist.append(vlist[i])
    #Find the vertices that are visible
    visible_vertices = [sortedvlist[0][4],sortedvlist[1][4]\
                        ,sortedvlist[2][4]] #List of int
    tflist = []
    for v in sortedvlist:
        tflist.append(TurtleFit(v))
    poly = [tflist[0],tflist[1],tflist[2]]
    for i in range(3,8):
        if in_poly(tflist[i], poly):
            break
        else:
            visible_vertices.append(sortedvlist[i][4])
            poly.append(tflist[i])
            poly = convex_hull(poly)
    #return the list of faces that are visible
    returnlist = []
    if 0 in visible_vertices and 1 in visible_vertices:
        if 2 in visible_vertices and 3 in visible_vertices:
            returnlist.append('xn')
        if 4 in visible_vertices and 5 in visible_vertices:
            returnlist.append('zp')
    if 4 in visible_vertices and 7 in visible_vertices:
        if 0 in visible_vertices and 3 in visible_vertices:
            returnlist.append('yp')
        if 5 in visible_vertices and 6 in visible_vertices:
            returnlist.append('xp')
    if 2 in visible_vertices and 6 in visible_vertices:
        if 1 in visible_vertices and 5 in visible_vertices:
            returnlist.append('yn')
        if 3 in visible_vertices and 7 in visible_vertices:
            returnlist.append('zn')
    return returnlist
        
def convex_hull(plist):
    """(list of (x,y)) --> list of (x,y)

    Find the polygon that encompasses the points specified.

    >>> convex_hull([(0,10),(10,0),(10,10),(0,0)])
    [(0, 0), (0, 10), (10, 10), (10, 0)]

    >>> convex_hull([(0,10),(10,0),(10,10),(0,0),(5,5)])
    [(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)]

    >>> convex_hull([(0,10),(10,0),(10,10),(0,0),(5,5),(5,15),(15,5),(15,15)])
    [(0, 0), (0, 10), (5, 15), (15, 15), (15, 5), (10, 0), (0, 0)]
    """
    length = len(plist)
    #find point with smallest x and y
    xlist = []
    ylist = []
    for c in plist:
        xlist.append(c[0])
        ylist.append(c[1])
    xmin = min(xlist)
    ymin = min(ylist)
##    pointOnHull = leftmost and lowest point in S
    returnlist = []
    pointOnHull = (xmin,ymin)
    firstpointOnHull = pointOnHull
    length-=1
    returnlist.append(pointOnHull)
    X = pointOnHull[0]
    Y = pointOnHull[1]
    endpoint = plist[1]
    for j in range(0,length):
        if endpoint == pointOnHull or \
           atan2((plist[j][1]-Y),(plist[j][0]-X))> atan2((endpoint[1]-Y),(endpoint[0]-X)):
            endpoint = plist[j]
    returnlist.append(endpoint)
    plist.remove(endpoint)
    length-=1
    pointOnHull = endpoint
    i = 1
##    repeat
##      P[i] = pointOnHull
    while returnlist[i] != firstpointOnHull and length >0:
        X = pointOnHull[0]
        Y = pointOnHull[1]
        endpoint = plist[0]
        #initial endpoint for a candidate edge on the hull
        for j in range(0,length+1):
            newangle = atan2((plist[j][1]-Y),(plist[j][0]-X))
            oldangle = atan2((endpoint[1]-Y),(endpoint[0]-X))
            if endpoint == pointOnHull or newangle > oldangle:
                #found greater left turn, update endpoint
                endpoint = plist[j]
        newY = endpoint[1]
        if newY < Y:
            break
        i += 1
        returnlist.append(endpoint)
        plist.remove(endpoint)
        length-=1
        pointOnHull = endpoint
        
    #second while loop handles second half of loop (when the Y starts decreasing)
    while returnlist[i] != firstpointOnHull and length >0:
        X = pointOnHull[0]
        Y = pointOnHull[1]
        endpoint = plist[0]
        #initial endpoint for a candidate edge on the hull
        for j in range(0,length+1):
            newangle = atan2fixed((plist[j][1]-Y),(plist[j][0]-X))
            oldangle = atan2fixed((endpoint[1]-Y),(endpoint[0]-X))
            if endpoint == pointOnHull or newangle > oldangle:
                #found greater left turn, update endpoint
                endpoint = plist[j]
        i += 1
        returnlist.append(endpoint)
        plist.remove(endpoint)
        length-=1
        pointOnHull = endpoint
    return returnlist

def atan2fixed(y,x):
    """(float,float) --> float

    Same as atan() in math module, except always returns positive angle value.

    >>> atan2fixed(1,1)
    0.7853981633974483

    >>> atan2fixed(-1,1)
    5.497787143782138

    Note: answer in radians
    """
    angle = atan2(y,x)
    if angle >= 0:
        return angle
    else:
        return 2*pi + angle
    
def in_poly(xytup, poly):
    """((x,y),list of (x,y)) --> bool

    reutrn whether or not the given point is in the polygon specified by the
    list of coordinates.
    """
    x = xytup[0]
    y = xytup[1]
    
    length = len(poly)
    i = 0
    j = length - 1
    c = False
    for i in range(length):
            if  ((poly[i][1] > y) != (poly[j][1] > y)) and \
                    (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
                c = not c
            j = i
    return c
    
def org_to_rot(space):
    """(space) -> Nonetype

    Converts original Space to roatated space. Allows new objects to be oriented
    with the screen.
    """
    global trans_list_one

    for trans in trans_list_one:
        print(trans)
        if trans[0] == 'r':
            space.rotate(trans[1][0],trans[1][1],trans[1][2])
        if trans[0] == 'm':
            space.translate(trans[1][0],trans[1][1],trans[1][2])
        if trans[0] == 'z':
            space.dilate(trans[1])
    
def main():
    
    def refreshHandler():
        global rot_list
        global trans_list
        global obj_list
        global oselected
        global sodict
        global axis_dic
        global trans_list_one
        forward_list = axis_display(axis_dic['vertices'])
        back_list = []
        t.clear()
        for s in axis_dic:
            if s != 'vertices' and s not in forward_list:
                back_list.append(s)
        t.pensize(1)
        for s in back_list:
            axis_dic[s].graph(screen,t)
        t.pensize(20)
        for space in rot_list:
            space.graph(screen,t)
        t.pensize(1)
        for s in forward_list:
            axis_dic[s].graph(screen,t)
        print('DONE')
        screen.update()
    def undoHandler():
        None 
    def loadHandler():
        None
    def saveHandler():
        None
################################################################################
    def rupHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('rup')
        for s in rot_list:
            s.rotate_x(-20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_x(-20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def rdownHandler():
        global rot_list
        global trans_list
        global axis_dic
        global trans_list_one
        trans_list.append('rdown')
        for s in rot_list:
            s.rotate_x(20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_x(20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()        
    def rleftHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('rleft')
        for space in rot_list:
            space.rotate_z(-20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_z(-20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def rrightHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('rright')
        for space in rot_list:
            space.rotate_z(20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_z(20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def rtupHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('rtup')
        for space in rot_list:
            space.rotate_y(-20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_y(-20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def rtdownHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('rtdown')
        for space in rot_list:
            space.rotate_y(20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].rotate_y(20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def zoutHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('zout')
        for space in rot_list:
            space.dilate(0.75)
            t.clear()
        for s in axis_dic:
            axis_dic[s].dilate(0.75)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def zinHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('zin')
        for space in rot_list:
            space.dilate(1.25)
            t.clear()
        for s in axis_dic:
            axis_dic[s].dilate(1.25)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def mupHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('mup')
        for space in rot_list:
            space.translate(0,0,20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(0,0,20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def mdownHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('mdown')
        for space in rot_list:
            space.translate(0,0,-20)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(0,0,-20)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def mleftHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('mleft')
        for space in rot_list:
            space.translate(-20,0,0)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(-20,0,0)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def mrightHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('mright')
        for space in rot_list:
            space.translate(20,0,0)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(20,0,0)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def minHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('min')
        for space in rot_list:
            space.translate(0,20,0)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(0,20,0)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
    def moutHandler():
        global rot_list
        global trans_list
        global obj_list
        global axis_dic
        global trans_list_one
        trans_list.append('mout')
        for space in rot_list:
            space.translate(0,-20,0)
            t.clear()
        for s in axis_dic:
            axis_dic[s].translate(0,-20,0)
            t.clear()
        trans_list_update(trans_list,trans_list_one)
        refreshHandler()
#################################################################################
    def moxpHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(40,0,0)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def moxnHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(-40,0,0)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def moypHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(0,40,0)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def moynHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(0,-40,0)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def mozpHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(0,0,40)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def moznHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        obj_list[oselected].translate(0,0,-40)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def reflectHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        Xprime = XprimeVar.get()
        Yprime = YprimeVar.get()
        Zprime = XprimeVar.get()
        if Xprime =='' or Yprime ==''or Zprime =='':
            obj_list[oselected].reflect(rfVar.get())
        else:
            obj_list[oselected].reflect(rfVar.get(),int(Xprime),int(Yprime),int(Zprime))
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def rotateHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        Xprime = XprimeVar.get()
        Yprime = YprimeVar.get()
        Zprime = XprimeVar.get()
        if Xprime =='' or Yprime ==''or Zprime =='':
            obj_list[oselected].rotate(int(rXVar.get()),int(rYVar.get()),int(rZVar.get()))
        else:
            obj_list[oselected].rotate(int(rXVar.get()),int(rYVar.get()),int(rZVar.get()),int(Xprime),int(Yprime),int(Zprime))
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def dilateHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        Xprime = XprimeVar.get()
        Yprime = YprimeVar.get()
        Zprime = XprimeVar.get()
        if Xprime =='' or Yprime ==''or Zprime =='':
            obj_list[oselected].dilate(float(dVar.get()))
        else:
            obj_list[oselected].dilate(float(dVar.get()),int(Xprime),int(Yprime),int(Zprime))
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def csHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        Xprime = XprimeVar.get()
        Yprime = YprimeVar.get()
        Zprime = XprimeVar.get()
        if Xprime =='' or Yprime ==''or Zprime =='':
            obj_list[oselected].compress(float(csXVar.get()),float(csYVar.get()),float(csZVar.get()))
        else:
            obj_list[oselected].compress(float(csXVar.get()),float(csYVar.get()),float(csZVar.get()),int(Xprime),int(Yprime),int(Zprime))
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def moHandler():
        global rot_list
        global obj_list
        global oselected
        global sodict
        mX = int(mXVar.get())
        mY = int(mYVar.get())
        mZ = int(mZVar.get())
        obj_list[oselected].translate(mX,mY,mZ)
        s = obj_list[oselected]
        sone = deepcopy(s)
        org_to_rot(sone)
        rot_list[oselected] = sone
        refreshHandler()
    def aoHandler():
        global obj_list
        global rot_list
        Z1 = Z1Var.get()
        increment = 20
        t.penup()
        Max = int(MaxVar.get())
        Min = int(MinVar.get())
        matrix = GraphList(Z1, Max,Min,increment)
        matrixone = deepcopy(matrix)
        rot_list.append(matrix)
        obj_list.append(matrixone)
        #select object
        oselected = len(rot_list)-1
        refreshHandler()
    def soHandler():
        global oselected
        oselected = sodict[soVar.get()]

    def clearHandler():
        global obj_list
        global rot_list
        global axis
        global trans_list
        global obj_list
        global sodict
        global axis_dic
        global trans_list_one
        t.clear()
        obj_list = []
        rot_list = []
        trans_list = []
        trans_list_one = []
        sodict = {}
        axis_dic_ini()
        
    root = tkinter.Tk()

    cv = tkinter.Canvas(root, width = 600, height=600)
    cv.pack(side=tkinter.LEFT)
    root.title("Graph")

    t = turtle.RawTurtle(cv)
    screen = t.getscreen()
    screen.setworldcoordinates(0,0,600,600)

    frame = tkinter.Frame(root)
    frame.pack(side = tkinter.RIGHT, fill = tkinter.BOTH)

    moveFrame = tkinter.Frame(frame,relief=tkinter.RIDGE, borderwidth=2)
    moveFrame.pack(side = tkinter.TOP,fill = tkinter.BOTH)

    rotateLabel = tkinter.Label(moveFrame, text = "Rotate")
    rotateLabel.pack(side=tkinter.TOP)

    rIFrame = tkinter.Frame(moveFrame)
    rIFrame.pack(side = tkinter.TOP)
    rIIFrame = tkinter.Frame(moveFrame)
    rIIFrame.pack(side = tkinter.TOP)
    rIIIFrame = tkinter.Frame(moveFrame)
    rIIIFrame.pack(side = tkinter.TOP)

    rupButton = tkinter.Button(rIFrame, text='^', command=rupHandler)
    rupButton.pack(side =tkinter.LEFT)
    rleftButton = tkinter.Button(rIIFrame, text="<", command=rleftHandler)
    rleftButton.pack(side =tkinter.LEFT)
    rrightButton = tkinter.Button(rIIFrame, text=">", command=rrightHandler)
    rrightButton.pack(side =tkinter.LEFT)
    rtdownButton = tkinter.Button(rIIFrame, text="v", command=rtdownHandler)
    rtdownButton.pack(side =tkinter.LEFT)
    rtupButton = tkinter.Button(rIIFrame, text='^', command=rtupHandler)
    rtupButton.pack(side =tkinter.LEFT)
    rdownButton = tkinter.Button(rIIIFrame, text="v", command=rdownHandler)
    rdownButton.pack(side =tkinter.LEFT)

    zFrame = tkinter.Frame(moveFrame)
    zFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    zoomLabel = tkinter.Label(zFrame, text = "Zoom")
    zoomLabel.pack(side=tkinter.LEFT)
    zinButton = tkinter.Button(zFrame, text="+", command=zinHandler)
    zinButton.pack(side =tkinter.LEFT)
    zoutButton = tkinter.Button(zFrame, text='-', command=zoutHandler)
    zoutButton.pack(side =tkinter.LEFT)

    moveLabel = tkinter.Label(moveFrame, text = "Move")
    moveLabel.pack(side=tkinter.TOP)

    mIFrame = tkinter.Frame(moveFrame)
    mIFrame.pack(side = tkinter.TOP)
    mIIFrame = tkinter.Frame(moveFrame)
    mIIFrame.pack(side = tkinter.TOP)
    mIIIFrame = tkinter.Frame(moveFrame)
    mIIIFrame.pack(side = tkinter.TOP)

    mupButton = tkinter.Button(mIFrame, text='^', command=mupHandler)
    mupButton.pack(side =tkinter.LEFT)
    mleftButton = tkinter.Button(mIIFrame, text="<", command=mleftHandler)
    mleftButton.pack(side =tkinter.LEFT)
    mrightButton = tkinter.Button(mIIFrame, text=">", command=mrightHandler)
    mrightButton.pack(side =tkinter.LEFT)
    moutButton = tkinter.Button(mIIFrame, text="v", command=moutHandler)
    moutButton.pack(side =tkinter.LEFT)
    minButton = tkinter.Button(mIIFrame, text='^', command=minHandler)
    minButton.pack(side =tkinter.LEFT)
    mdownButton = tkinter.Button(mIIIFrame, text="v", command=mdownHandler)
    mdownButton.pack(side =tkinter.LEFT)

    ooFrame = tkinter.Frame(frame,relief=tkinter.RIDGE, borderwidth=2)
    ooFrame.pack(side=tkinter.TOP)
    ooLabel = tkinter.Label(ooFrame, text = 'Object Options')
    ooLabel.pack(side=tkinter.TOP)
    
    centerFrame = tkinter.Frame(ooFrame)
    centerFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    centerLabel = tkinter.Label(centerFrame, text = 'Center')
    centerLabel.pack(side=tkinter.LEFT)
    XprimeLabel = tkinter.Label(centerFrame, text="X'")
    XprimeLabel.pack(side=tkinter.LEFT)
    XprimeVar = tkinter.StringVar()
    XprimeEntry = tkinter.Entry(centerFrame, textvariable=XprimeVar,width=5)
    XprimeEntry.pack(side=tkinter.LEFT)
    YprimeLabel = tkinter.Label(centerFrame, text="Y'")
    YprimeLabel.pack(side=tkinter.LEFT)
    YprimeVar = tkinter.StringVar()
    YprimeEntry = tkinter.Entry(centerFrame, textvariable=YprimeVar,width=5)
    YprimeEntry.pack(side=tkinter.LEFT)
    ZprimeLabel = tkinter.Label(centerFrame, text="Z'")
    ZprimeLabel.pack(side=tkinter.LEFT)
    ZprimeVar = tkinter.StringVar()
    ZprimeEntry = tkinter.Entry(centerFrame, textvariable=ZprimeVar,width=5)
    ZprimeEntry.pack(side=tkinter.LEFT)

    rotateFrame = tkinter.Frame(ooFrame)
    rotateFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    rotateLabel = tkinter.Label(rotateFrame, text = 'Rotate')
    rotateLabel.pack(side=tkinter.LEFT)
    rXLabel = tkinter.Label(rotateFrame, text="X")
    rXLabel.pack(side=tkinter.LEFT)
    rXVar = tkinter.StringVar()
    rXEntry = tkinter.Entry(rotateFrame, textvariable=rXVar,width=5)
    rXEntry.pack(side=tkinter.LEFT)
    rYLabel = tkinter.Label(rotateFrame, text="Y")
    rYLabel.pack(side=tkinter.LEFT)
    rYVar = tkinter.StringVar()
    rYEntry = tkinter.Entry(rotateFrame, textvariable=rYVar,width=5)
    rYEntry.pack(side=tkinter.LEFT)
    rZLabel = tkinter.Label(rotateFrame, text="Z")
    rZLabel.pack(side=tkinter.LEFT)
    rZVar = tkinter.StringVar()
    rZEntry = tkinter.Entry(rotateFrame, textvariable=rZVar,width=5)
    rZEntry.pack(side=tkinter.LEFT)
    rgoButton = tkinter.Button(rotateFrame, text='Go', command=rotateHandler)
    rgoButton.pack(side=tkinter.LEFT)

    reflectFrame = tkinter.Frame(ooFrame)
    reflectFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    reflectLabel = tkinter.Label(reflectFrame, text = 'Reflect')
    reflectLabel.pack(side=tkinter.LEFT)
    rfVar = tkinter.StringVar()
    Rxy = tkinter.Radiobutton(reflectFrame, text="xy", variable=rfVar, value='xy')
    Rxy.pack(side=tkinter.LEFT)
    Rxz = tkinter.Radiobutton(reflectFrame, text="xz", variable=rfVar, value='xz')
    Rxz.pack(side=tkinter.LEFT)
    Ryz = tkinter.Radiobutton(reflectFrame, text="yz", variable=rfVar, value='yz')
    Ryz.pack(side=tkinter.LEFT)
    rfgoButton = tkinter.Button(reflectFrame, text='Go', command=reflectHandler)
    rfgoButton.pack(side=tkinter.LEFT)

    dilateFrame = tkinter.Frame(ooFrame)
    dilateFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    dilateLabel = tkinter.Label(dilateFrame, text = 'Dilate   Factor')
    dilateLabel.pack(side=tkinter.LEFT)
    dVar = tkinter.StringVar()
    dEntry = tkinter.Entry(dilateFrame, textvariable=dVar,width=5)
    dEntry.pack(side=tkinter.LEFT)
    dgoButton = tkinter.Button(dilateFrame, text='Go', command=dilateHandler)
    dgoButton.pack(side=tkinter.LEFT)

    csFrame = tkinter.Frame(ooFrame)
    csFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    csLabel = tkinter.Label(csFrame, text = 'Compress/Stretch')
    csLabel.pack(side=tkinter.LEFT)
    csXLabel = tkinter.Label(csFrame, text="X")
    csXLabel.pack(side=tkinter.LEFT)
    csXVar = tkinter.StringVar()
    csXEntry = tkinter.Entry(csFrame, textvariable=csXVar,width=5)
    csXEntry.pack(side=tkinter.LEFT)
    csYLabel = tkinter.Label(csFrame, text="Y")
    csYLabel.pack(side=tkinter.LEFT)
    csYVar = tkinter.StringVar()
    csYEntry = tkinter.Entry(csFrame, textvariable=csYVar,width=5)
    csYEntry.pack(side=tkinter.LEFT)
    csZLabel = tkinter.Label(csFrame, text="Z")
    csZLabel.pack(side=tkinter.LEFT)
    csZVar = tkinter.StringVar()
    csZEntry = tkinter.Entry(csFrame, textvariable=csZVar,width=5)
    csZEntry.pack(side=tkinter.LEFT)
    csgoButton = tkinter.Button(csFrame, text='Go', command=csHandler)
    csgoButton.pack(side=tkinter.LEFT)

    moveFrame = tkinter.Frame(ooFrame)
    moveFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    moIFrame = tkinter.Frame(moveFrame)
    moIFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    moveLabel = tkinter.Label(moIFrame, text = 'Move')
    moveLabel.pack(side=tkinter.LEFT)
    mXLabel = tkinter.Label(moIFrame, text="X")
    mXLabel.pack(side=tkinter.LEFT)
    mXVar = tkinter.StringVar()
    mXEntry = tkinter.Entry(moIFrame, textvariable=mXVar,width=5)
    mXEntry.pack(side=tkinter.LEFT)
    mYLabel = tkinter.Label(moIFrame, text="Y")
    mYLabel.pack(side=tkinter.LEFT)
    mYVar = tkinter.StringVar()
    mYEntry = tkinter.Entry(moIFrame, textvariable=mYVar,width=5)
    mYEntry.pack(side=tkinter.LEFT)
    mZLabel = tkinter.Label(moIFrame, text="Z")
    mZLabel.pack(side=tkinter.LEFT)
    mZVar = tkinter.StringVar()
    mZEntry = tkinter.Entry(moIFrame, textvariable=mZVar,width=5)
    mZEntry.pack(side=tkinter.LEFT)
    mogoButton = tkinter.Button(moIFrame, text='Go', command=moHandler)
    mogoButton.pack(side=tkinter.LEFT)

    moFrame = tkinter.Frame(moveFrame)
    moFrame.pack(side=tkinter.TOP,fill = tkinter.BOTH)
    moxFrame = tkinter.Frame(moFrame)
    moxFrame.pack(side = tkinter.LEFT)
    moyFrame = tkinter.Frame(moFrame)
    moyFrame.pack(side = tkinter.LEFT)
    mozFrame = tkinter.Frame(moFrame)
    mozFrame.pack(side = tkinter.LEFT)
    moxLabel = tkinter.Label(moxFrame,text= 'X')
    moxLabel.pack(side =tkinter.LEFT)
    moyLabel = tkinter.Label(moyFrame,text= 'Y')
    moyLabel.pack(side =tkinter.LEFT)
    mozLabel = tkinter.Label(mozFrame,text= 'Z')
    mozLabel.pack(side =tkinter.LEFT)
    moxpButton = tkinter.Button(moxFrame, text='+', command=moxpHandler)
    moxpButton.pack(side =tkinter.LEFT)
    moxnButton = tkinter.Button(moxFrame, text="-", command=moxnHandler)
    moxnButton.pack(side =tkinter.LEFT)
    moypButton = tkinter.Button(moyFrame, text="+", command=moypHandler)
    moypButton.pack(side =tkinter.LEFT)
    moynButton = tkinter.Button(moyFrame, text="-", command=moynHandler)
    moynButton.pack(side =tkinter.LEFT)
    mozpButton = tkinter.Button(mozFrame, text='+', command=mozpHandler)
    mozpButton.pack(side =tkinter.LEFT)
    moznButton = tkinter.Button(mozFrame, text="-", command=moznHandler)
    moznButton.pack(side =tkinter.LEFT)

    aoFrame = tkinter.Frame(frame)
    aoFrame.pack(side=tkinter.TOP)
    aoVar = tkinter.StringVar()
    aochoices = ['Custom']
    aoOption = tkinter.OptionMenu(aoFrame, aoVar, *aochoices)
    aoOption.pack(side=tkinter.LEFT)
    aoButton = tkinter.Button(aoFrame,text="Add Object",command=aoHandler)
    aoButton.pack(side =tkinter.LEFT)
    
    
    screen.tracer(0)    
    
    Z1Lab = tkinter.Label(frame, text="Z1= ")
    Z1Lab.pack()

    Z1Var = tkinter.StringVar()
    Z1Entry = tkinter.Entry(frame, textvariable=Z1Var)
    Z1Entry.pack()

    MaxMinFrame = tkinter.Frame(frame)
    MaxMinFrame.pack(side=tkinter.TOP)
    maxLabel = tkinter.Label(MaxMinFrame,text= 'Max')
    maxLabel.pack(side =tkinter.LEFT)
    MaxVar = tkinter.StringVar()
    maxEntry = tkinter.Entry(MaxMinFrame, textvariable=MaxVar, width=5)
    maxEntry.pack(side=tkinter.LEFT)
    minLabel = tkinter.Label(MaxMinFrame,text= 'Min')
    minLabel.pack(side =tkinter.LEFT)
    MinVar = tkinter.StringVar()
    minEntry = tkinter.Entry(MaxMinFrame, textvariable=MinVar, width=5)
    minEntry.pack(side=tkinter.LEFT)

    clearButton = tkinter.Button(frame, text="Clear", command=clearHandler)
    clearButton.pack()
    root.mainloop()

main()
