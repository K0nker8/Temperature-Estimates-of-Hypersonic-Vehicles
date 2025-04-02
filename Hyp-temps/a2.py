import matplotlib.pyplot as plt
import numpy as np
import random
import math
import csv

def nearest_edge_square(x, y, length):
    '''Function to find the distance to and edge number of the closest edge from a point within a sqaure.

    Parameters:
    x (float): The x position of a point within the square
    y (float): The y position of a point within the square
    length (float): The length of the square
    
    Returns:
    square_edge (tuple(int, float)): A tuple containing the edge index and distance to the closest edge.
    '''
    #Defines the distances to each edge based on a square domain starting at 0.
    #(e.g: y represents the perpendular distance to y-axis origin, and hence the first edge).
    length_0 = y
    length_1 = length - x
    length_2 = length - y
    length_3 = x

    # Finds the minimum distance to an edge, and notes its index by first putting each into a list. 
    distances = [length_0,length_1,length_2,length_3]
    edge = distances.index(min(distances))
    dist = min(distances)
    square_edge = (edge, dist)
    return square_edge

def temperature_at_nearest_edge(x, y, length, edge, cnr_temps):
    ''' Function to find the temperature at the nearest edge of the square to the point,
    accounting for perpedicular distance of the point along the edge.
    
    Parameters:
    x (float): The x position of a point within the square.
    y (float): The y position of a point within the square.
    length (float): The length of the square.
    edge (int): The edge index closest to the point.
    cnr_temps(dict): Dictionary of temperatures at the 4 corners of the square.

    Returns:
    T_edge (float): The temperature at the edge closest to the point, accounting for perpendicular distance.
    
    '''
    #Finds the fractional distance along the closest line based on the position in the other dimension.
    if edge == 0 or edge == 2:
        fract_dist = x / length
    elif edge == 1 or edge == 3:
        fract_dist = y / length
    #Weights the temeprature at each corner by the fractional distance of the point along the line.
    if edge == 0:
        T_edge = (1 - fract_dist) * cnr_temps["sw"] + fract_dist * cnr_temps["se"]
        
    elif edge == 1:
        T_edge = (1 - fract_dist) * cnr_temps["se"] + fract_dist * cnr_temps["ne"]

    elif edge == 2:
        T_edge = (1 - fract_dist) * cnr_temps["nw"] + fract_dist * cnr_temps["ne"]

    elif edge == 3:
        T_edge = (1 - fract_dist) * cnr_temps["sw"] + fract_dist * cnr_temps["nw"]
    return T_edge

def jump_on_circle(x, y, radius):
    '''
    Function to radomly move the point by a set radius at a random angle.
    Parameters:
    x (float): The x position of a point within the square.
    y (float): The y position of a point within the square.
    radius(float): The radius that the point is to move.

    Returns:
    x_circ(float): The new x cooridnate after moving a set distance
    y_circ(float): The new y cooridnate after moving a set distance   
    '''
    # Randomises a direction to move and then moves the point a distance radius.
    angle = random.uniform(-math.pi,math.pi)
    x_circ = x + radius * math.cos(angle)
    y_circ = y + radius * math.sin(angle)
    return (x_circ, y_circ)
    
def wos_temperature_estimate(x0, y0, length, cnr_temps, epsilon, number_walks):
    '''
    Function to get a walk on spheres temeprature estimate of temperature on a square domain.
    
    Parameters:
    x0 (float): The initial x posotion of a point before performing each walk.
    y0 (float): The initial y posotion of a point before performing each walk.
    length (float): The length of the square domain.
    cnr_temps (dict): Dictionary of temperatures at the 4 corners of the square.
    epsilon (float): Distance considered close enough to stop each walk.
    number_walks (int): The number of walks to be performed.
    
    Returns:
    average (float): The average temperature from the sum of all walks performed.
    numpy_temps (numpy array): Each temeprature value computed from each walk.
    
    '''
    x = x0
    y = y0
    temp_vals = []
    #Performes jump on circle until the point is with epsilon distance from the edge.
    for i in range(number_walks):
        while epsilon <= nearest_edge_square(x, y, length)[1]:
            edge_dist = nearest_edge_square(x, y, length)
            r = edge_dist[1]
            new_point = jump_on_circle(x, y, r)
            x = new_point[0]
            y = new_point[1]
        #Finds the temeprature at the closest edge to the new point and appends it to a list.
        #Resets to the orignal point at the end of the walk
        edge_dist = nearest_edge_square(x, y, length)
        new_temp = temperature_at_nearest_edge(x, y, length, edge_dist[0], cnr_temps)
        temp_vals.append(new_temp)
        x = x0
        y = y0
    temp_sum = 0
    #Finds an average temperature from the list of walks.
    for i in temp_vals:
        temp_sum += i
    average = temp_sum / len(temp_vals)
    numpy_temps = np.array(temp_vals)
    return (average, numpy_temps)

def plot_temperature_histogram(Ts, bins, range):
    '''
    Function to plot a histogram of the individual temepratures from each walk,
    grouped into counts based on temperature ranges.

    Parameters:
    Ts (numpy array): Temperature values from each walk
    bins (int): The number of bins to be used in the histogram
    range (tuple(float, float)): Tuple of the minimum and maximum temperature bounds for the histogram.

    Returns:
    None
    
    '''
    plt.hist(Ts, bins=bins, range=range)
    plt.xlabel('temperature, K')
    plt.ylabel('count')

class Vec2D():
    '''
    An object representing a two dimensional vector on a cartesian plane.
    '''

    def __init__(self, x, y):
        '''
        Function to initialise x and y psotions of the vector
        
        Parameters:
        self: The vector object
        x (float): The x cooridinate of the vector
        y (float): The y cooridinate of the vector

        Returns:
        None
        '''
        
        self.x = x
        self.y = y
        
    def __str__(self):
        '''
        Method to print a string of the x and y coordinates of the vector.
        Parameters:
        self: The vector object.
        Returns:
        (self.x, self.y) (tuple(str, str)): Strings of the x and y coordinates of the vector.
        '''
        return str((self.x,self.y))
    
    def __repr__(self):
        '''
        Method to print a string of the x and y coordinates of the vector showing that it is a Vec2D object.
        
        Parameters:
        self: The vector object.
        
        Returns:
        Vec2D(self.x, self.y) (str + tuple (str, str)): Strings of the x and y coordinates of the vector,
        showing that is is a Vec2D object.
        '''
        return "Vec2D" + str(Vec2D(self.x,self.y))
    
    def get_x(self):
        '''
        Gets the x coordiate of the vector.
        
        Parameters:
        self: The vector object.
        
        Returns:
        self.x (str): The x coordinate of the vector. 
        '''
        return self.x
    
    def get_y(self):
        '''
        Gets the y coordiate of the vector.
        
        Parameters:
        self: The vector object.
        
        Returns:
        self.y (str): The y coordinate of the vector. 
        '''
        return self.y
        
    def set_x(self, x):
        '''
        Sets the x coordiate of the vector to a new value.
        
        Parameters:
        self: The vector object.
        x (float): New x value to set to
        
        Returns:
        None
        '''
        self.x = x
        
    def set_y(self, y):
        '''
        Sets the y coordiate of the vector to a new value.
        
        Parameters:
        self: The vector object.
        y (float): New x value to set to
        
        Returns:
        None
        '''
        self.y = y
        
    def __add__(self, b):
        '''
        Method to add the self vector to another vector.
        
        Parameters:
        self: The vector object.
        b (Vec2D(float, float)): Vec2D object of the vector to add to.
        
        Returns:
        self.x + b.x (Vec2D(float)): Sum of the x components of the two vectors.
        self.y + b.y (Vec2D(float)): Sum of the y components of the two vectors.
        
        '''
        return Vec2D(self.x + b.x, self.y + b.y)
        
    def __sub__(self, b):
        '''
        Method to substract another vector from the self vector. 
        Parameters:
        self: The vector object.
        b (Vec2D(float, float)): Vec2D object of the vector being subtracted.
        
        Returns:
        self.x - b.x (Vec2D(float)): X component of the b vector subtracted from the self vector.
        self.y - b.y (Vec2D(float)): X component of the a vector subtracted from the self vector.
        '''
        return Vec2D(self.x - b.x, self.y - b.y)
    
    def __mul__(self, scalar):
        
        '''
        Method to multiply the self vector to another vector on the left.
        
        Parameters:
        self: The vector object.
        scalar (float): A scalar multiple 
        
        Returns:
        self.x * scalar (Vec2D(float)): Multiplication of the x components of the two vectors on the left.
        self.y * scalar (Vec2D(float)): Multiplication of the y components of the two vectors on the left.
        
        '''
        return Vec2D(self.x * scalar, self.y * scalar)
    def __rmul__(self, scalar):
        '''
        Method to multiply the self vector to another vector on the right.
        
        Parameters:
        self: The vector object.
        scalar (float): A scalar multiple 
        
        Returns:
        scalar * self.x  (Vec2D(float)): Multiplication of the x components of the two vectors on the right.
        scalar * self.y (Vec2D(float)): Multiplication of the y components of the two vectors on the right.
        
        '''
        return Vec2D(scalar * self.x, scalar * self.y)
        
    
    def length(self):
        '''
        Method to find the length of a vector

        Parameters:
        self: The vector object

        Returns:
        sqrt(self.x * self.x + self.y * self.y) (float): The square root of the square of
        the x and y components of the vector.
        '''
        return math.sqrt(self.x * self.x + self.y * self.y)

    def dot(self, b):
        '''
        Method to find the dot product of the self vector and another vector.

        Parameters:
        self: The vector object
        b: (Vec2D(float, float)): Vec2D object of the vector being dotted with.

        Returns:
        self.x * b.x + self.y * b.y (float): The dot product of the self vector with the b vector
        '''
        return self.x * b.x + self.y * b.y
    
def closest_point_on_line_segment(p, a, b):
    '''
    Method to find the closest point to a vector p on a line segment formed by the vectors a and b.

    Parameters:
    p (Vec2D(float, float)): Vec2D coordinates of the point within the object
    a:(Vec2D(float, float)): Vec2D coordinates of the first point on the line
    b:(Vec2D(float, float)): Vec2D coordinates of the second point on the line

    Returns:
    c: (Vec2D(float, float)): Vec2D coordinates of the closest point to p formed by the vectors a and b.
    
    '''
    #Finds the distance bertween the start and end point of the edge as u
    #Finds the projection of the vector from a to p onto u
    #Finds a point c, which is the linear weighting using t bertween a and b. 
    u = b - a
    t = (p - a).dot(u) / (u.dot(u))
    if t < 0:
        t = 0.0
    elif t > 1:
        t = 1.0
    c = ((1 - t) * a) + (t * b)
    return c

class Boundary():
    '''
    A polygon object with defined edges in which the jump on circles method will be performed.
    '''
    def __init__(self, filename):
        '''
        Initialising method that opens a data file and stores its
        x, y and repective temperature values as the start and endpoints of each edge.
        
        Parameters:
        self: The polygon object
        filename(.dat file): The file containing the coordinates and temperatures of each point in the polygon

        Returns:
        None
        '''
        #Splits each line of the file into a list of three values, and then appends each value into ist own list
        with open(filename, 'r') as file:
            self.x_vals = []
            self.y_vals = []
            self.temp_vals = []
            for line in file:
                split_line = line.split()
                self.x_vals.append(float(split_line[0]))
                self.y_vals.append(float(split_line[1]))
                self.temp_vals.append(float(split_line[2]))
                            
    def distance_to_closest_edge(self, p):
        '''
        Method to find the distance to the closest edge given a point within the polygon and the corresponding edge index.
        
        Parameters:
        self: the polygon object
        p (Vec2D(float, float)): Vec2D coordinates of the point within the polygon

        Returns:
        closest_distance_index (float): Index of the closest line segment to the p vector
        corresponding to the order listed in the file
        closest_distance (float): Distance bertween the p vector and the closest line segment.
        '''
        edge_distances = []
        #Iterates through each edge defined by the start and end points. Finds the closest point on each line segment to the point.
        #Finds the distance bertween the point and the closest point on that edge.
        for count in range(len(self.x_vals) - 1):
             segment_position = (closest_point_on_line_segment(p, Vec2D(self.x_vals[count], self.y_vals[count]), Vec2D(self.x_vals[count+1], self.y_vals[count+1])))
             edge_distance = (segment_position - p).length()
             edge_distances.append(edge_distance)
        #Iterates through each disatance and finds the shortest and its corresponding index.
        closest_distance = min(edge_distances)
        closest_distance_index = edge_distances.index(min(edge_distances))
        return (closest_distance_index, closest_distance)
             
        
    def temperature_at_edge(self, p, edge_index):
        '''
        Method to find the temperature at the closest point on the edge to the point.

        Parameters:
        self: the polygon object
        p (Vec2D(float, float)): Vec2D coordinates of the point within the polygon
        edge_index (int): Index of the edge returned from distance_to_closest_edge

        Returns:
        t_edge (float): Temperature at the closest point on the edge to the point.
        
        '''
        p_x = p.get_x()
        p_y = p.get_y()

        # The values at the start of the edge
        x1 = self.x_vals[edge_index]
        y1 = self.y_vals[edge_index]
        temp1 = self.temp_vals[edge_index]
        # The values at the end of the edge
        x2 = self.x_vals[edge_index + 1]
        y2 = self.y_vals[edge_index + 1]
        temp2 = self.temp_vals[edge_index + 1]

        a = Vec2D(x1, y1)
        b = Vec2D(x2, y2)
        p = Vec2D(p_x, p_y)
        
        #Finds the distance bertween the start and end point of the edge as u.
        #Finds the projection of the vector from a to p onto u.
        #Finds a point t_edge, which is the linear weighting using t bertween a and b. 
        u = b - a
        t = (p - a).dot(u) / (u.dot(u))
        t_edge = (1-t) * temp1 + t * temp2

        return t_edge

    def wos_temperature_estimate(self, p0, epsilon, number_walks):
        '''
        Method to run the walk on spheres multiple times from a set point within the polygon
        until it reaches a point less than epsilon. Finds the avearge tempearture from the walks performed.

        Parameters:
        self: the polygon object
        p0 (Vec2D(float, float)): Vec2D coordinates of the intial point within the polygon
        epsilon (float): Distance considered to be samll enough for the point to be considered at the line.
        number_walks (int): The number of walks to be performed

        Returns:
        average (float): Average of the temperatures from performing all walks.
        numpy_temps (numpy array): The resultant temperatures after performing each walk.
        
        '''
        temp_vals =  []
        p = p0
        #For a set number of iterations, performs the walk method until the point is within eplison distance of an edge.
        for i in range(number_walks):  
            while epsilon < self.distance_to_closest_edge(p)[1]:
                #Sets the distance to move as the distance to the current closest edge.
                radius = self.distance_to_closest_edge(p)[1]
                #Uses jump on circle to set a new position
                new_point = jump_on_circle(p.get_x(), p.get_y(), radius)
                p = Vec2D(new_point[0], new_point[1])
            #Find the temeprature at nearest edge to point after the walk is complete. 
            new_temp = self.temperature_at_edge(p, self.distance_to_closest_edge(p)[0])
            temp_vals.append(new_temp)
            p = p0
        temp_sum = 0
        #Finds a sum of the temepratures found from all walks and takes an average.
        for i in temp_vals:
            temp_sum += i
        average = temp_sum / len(temp_vals)
        numpy_temps = np.array(temp_vals)
        return (average, numpy_temps)
    
def plot_temperature_convergence(Ts):
    '''
    Function to plot a running average of the temperatures from each walk and the
    temeprature value for each with the walk number.
    
    Parameters:
    Ts: (numpy array): Array of resultant temepratures after each walk

    Returns:
    None
    '''
    count = 1
    sum = 0
    walk_numbers = []
    moving_averages = []
    #Calculates a running average for each temeprature value in Ts.
    for i in Ts:
        walk_numbers.append(count)
        sum += float(Ts[count - 1])
        new_average = sum / count
        moving_averages.append(new_average)
        count += 1

    #Setting up moving average.
    walk_number = np.array(walk_numbers)
    temp_estimates = np.array(Ts)
    moving_averages = np.array(moving_averages)

    #Plots moving average and other data.
    fig,ax = plt.subplots()
    ax.plot(walk_number, temp_estimates, marker='.', linestyle='None', color='blue', label="edge temperature per walk")
    ax.plot(walk_number, moving_averages, color='red', label="running average temperature")
    
    
    plt.xlabel("walk number")
    plt.ylabel("temperature, k")
    plt.grid()
    plt.legend()
    plt.show()
                 
            
        



        
        
        
            
        
        

            
           
       
        
        
    
        
        
        
        
        
    

    
    
    
    
        
        
        
        
        
       
    
    
    
            
