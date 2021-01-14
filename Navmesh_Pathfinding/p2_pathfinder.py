import math
from heapq import heappush, heappop

# MADE BY:
#   Herman Wu, hwwu
#   Justin Lao, jlao3

def find_path (source_point, destination_point, mesh):
    """
    Searches for a path from source_point to destination_point through the mesh

    Args:
        source_point: starting point of the pathfinder
        destination_point: the ultimate goal the pathfinder must reach
        mesh: pathway constraints the path adheres to

    Returns:

        A path (list of points) from source_point to destination_point if exists
        A list of boxes explored by the algorithm
    """
    src_box = None
    des_box = None

    for x1, x2, y1, y2 in mesh['boxes']:
        if x1 <= source_point[0] <= x2 and y1 <= source_point[1] <= y2:
            src_box = (x1, x2, y1, y2)
        if x1 <= destination_point[0] <= x2 and y1 <= destination_point[1] <= y2:
            des_box = (x1, x2, y1, y2)

    if src_box is None and des_box is None:
        return noPathFound(debugMessage = "src_box and des_box are not defined")
    elif src_box is None:
        return noPathFound(debugMessage = "src_box is not defined")
    elif des_box is None:
        return noPathFound(debugMessage = "des_box is not defined")

    # path, boxes = bfs( mesh, source_point, destination_point, src_box, des_box )
    # path, boxes = dijkstras( mesh, source_point, destination_point, src_box, des_box )
    # path, boxes = a_star( mesh, source_point, destination_point, src_box, des_box )
    path, boxes = bidirectional_a_star( mesh, source_point, destination_point, src_box, des_box )

    return path, boxes

# Returns euclidean distance between two points
def getDistance(a, b):
    return math.sqrt( (a[0] - b[0])**2 + (a[1] - b[1])**2 )

# Calculating x1, x2, y1, y2 since we're working with rectangles
def calculate_distance(curr_box_coords, destination_point, adj_box):
    x1, x2 = adj_box[0], adj_box[1]
    next_x = min(x2 - 1, max(x1, curr_box_coords[0]))

    y1, y2 = adj_box[2], adj_box[3]
    next_y = min(y2 - 1, max(y1, curr_box_coords[1]))

    point_in_neighbor = (next_x, next_y)

    total_dist = getDistance( point_in_neighbor, curr_box_coords )
    
    dist_to_goal = getDistance( point_in_neighbor, destination_point )

    return total_dist, dist_to_goal, point_in_neighbor

# Returns the midpoint of box
def boxMidpoint( box ):
    xMid = ( box[0] + box[1] ) / 2
    yMid = ( box[2] + box[3] ) / 2
    return ( xMid, yMid )

# Returns empty data structures
def noPathFound( boxes = {}, debugMessage = '' ):
    print("No path! " + debugMessage)
    return [], boxes

# Given a dictionary of paths for boxes (came_from), return a path of coordinates
def getCoordPathWithMidpoints( came_from, source_point, destination_point, src_box, des_box ):
    coord_path = [destination_point]            # Begin construction of the coord_path
    current_box = came_from[des_box]

    if current_box == None:                     # If src_box == des_box
        coord_path.insert(0, source_point )
    elif current_box == src_box:                # If src_box neighbors des_box
        # Likely will need to do some additional things here when we swap from midpoints
        coord_path.insert(0, source_point )
    else:
        while current_box is not src_box:
            coord_path.insert(0, boxMidpoint( current_box ) )
            current_box = came_from[current_box]
        coord_path.insert(0, source_point )
    
    return coord_path

def getCoordPathFromDict( came_from, box_coord, source_point, destination_point, src_box, des_box ):
    coord_path = [destination_point]            # Begin construction of the coord_path
    current_box = came_from[des_box]

    if current_box == None:                     # If src_box == des_box
        coord_path.insert(0, source_point )
    elif current_box == src_box:                # If src_box neighbors des_box
        coord_path.insert(0, box_coord[ des_box ] )
        coord_path.insert(0, box_coord[ current_box ] )
        coord_path.insert(0, source_point )
    else:
        coord_path.insert(0, box_coord[ des_box ] )
        while current_box is not src_box:
            coord_path.insert(0, box_coord[ current_box ] )
            current_box = came_from[current_box]
        coord_path.insert(0, source_point )
    
    return coord_path

def getCoordPathForBiDirAStar( forward_came_from, backward_came_from, forward_box_coord, backward_box_coord, source_point, destination_point, intersection_box, src_box, des_box ):
    forward_coord_path = []
    current_box = forward_came_from[intersection_box]

    if current_box == None:                     # If src_box == intersection_box
        forward_coord_path.insert(0, source_point )
    elif current_box == src_box:                # If src_box neighbors intersection_box
        forward_coord_path.insert(0, forward_box_coord[ intersection_box ] )
        forward_coord_path.insert(0, forward_box_coord[ current_box ] )
        forward_coord_path.insert(0, source_point )
    else:
        forward_coord_path.insert(0, forward_box_coord[ intersection_box ] )
        while current_box is not src_box:
            forward_coord_path.insert(0, forward_box_coord[ current_box ] )
            current_box = forward_came_from[current_box]
        forward_coord_path.insert(0, source_point )

    backward_coord_path = []            # Begin construction of the coord_path
    current_box = backward_came_from[intersection_box]

    if current_box == None:                     # If src_box == des_box
        backward_coord_path.insert(0, destination_point )
    elif current_box == src_box:                # If src_box neighbors des_box
        backward_coord_path.insert(0, backward_box_coord[ intersection_box ] )
        backward_coord_path.insert(0, backward_box_coord[ current_box ] )
        backward_coord_path.insert(0, destination_point )
    else:
        backward_coord_path.insert(0, backward_box_coord[ intersection_box ] )
        while current_box is not des_box:
            backward_coord_path.insert(0, backward_box_coord[ current_box ] )
            current_box = backward_came_from[current_box]
        backward_coord_path.insert(0, destination_point )
    
    backward_coord_path.reverse()
    
    return forward_coord_path + backward_coord_path

def bfs( mesh, source_point, destination_point, src_box, des_box ):
    boxes = {}          # Dict tracking all boxes visited
    priorityQ = []      # Heap tracking what box to investigate next

    # back pointer
    came_from = {src_box: None}

    # coordinates of each box pointer
    box_coord = {src_box: source_point}

    # keeps track of cost
    cost_so_far = {src_box: 0}

    heappush(priorityQ, (0, src_box))

    while priorityQ:
        distance, curr_box = heappop(priorityQ)

        # Finish the search if we reached the destination
        if curr_box == des_box:
            break

        for adj_box in mesh['adj'][curr_box]:
            boxes[adj_box] = curr_box       # making the pointer equal to curr box
            if adj_box not in came_from:
                heappush(priorityQ, (0, adj_box))
                came_from[adj_box] = curr_box

    # Checks if destination exists in came_from
    # If it doesn't, then we didn't reach it
    if des_box not in came_from.keys():
        return noPathFound(boxes, "BFS could not reach destination")

    path = getCoordPathWithMidpoints( came_from, source_point, destination_point, src_box, des_box )

    return  path, boxes

def dijkstras( mesh, source_point, destination_point, src_box, des_box ):
    path = []
    boxes = {}

    # back pointer
    came_from = {src_box: None}

    # coordinates of each box pointer
    box_coord = {src_box: source_point}

    # keeps track of cost
    cost_so_far = {src_box: 0}

    priorityQ = [(0, src_box)]

    while priorityQ:
        current_dist, current_box = heappop(priorityQ)

        # Finish the search if we reached the destination
        if current_box == des_box:
            break
        
        # Calculate cost from current box to all the adjacent ones
        for adj_box in mesh['adj'][current_box]:
            boxes[adj_box] = current_box

            # Finds edge range
            x_range = [ max(current_box[0], adj_box[0]), min(current_box[1], adj_box[1]) ]
            y_range = [ max(current_box[2], adj_box[2]), min(current_box[3], adj_box[3]) ]

            # Mid point of edge range
            mid_point = ( (x_range[0] + x_range[1]) * 0.5, (y_range[0] + y_range[1]) * 0.5 )

            """
            mid_slope = ( destination_point[1] - box_coord[current_box][1] ) / ( destination_point[0] - box_coord[current_box][0] )

            if max( x_range ) - min( x_range ) == 0:                                                                # If the border is vertical
                y = mid_slope * ( x_range[0] - box_coord[current_box][0] ) + box_coord[current_box][1]              # Find the y coordinate
            else:                                                                                                   # If the border is horizontal
                x = ( 1 / mid_slope ) * ( y_range[0] - box_coord[current_box][1] ) + box_coord[current_box][0]      # FInd the x coordinate
            """
            # points at the ends of the edge
            a_point = (x_range[0],y_range[0])
            b_point = (x_range[1],y_range[1])

            # Calculates cost from each edge point
            a_cost = getDistance( box_coord[current_box], a_point ) #+ getDistance( a_point, destination_point )
            b_cost = getDistance( box_coord[current_box], b_point ) #+ getDistance( b_point, destination_point )
            mid_cost = getDistance( box_coord[current_box], mid_point ) #+ getDistance( mid_point, destination_point )

            # Lowest cost is added to the queue
            if a_cost <= b_cost and a_cost <= mid_cost:
                edge_cost = a_cost
                edge_point = a_point
            elif b_cost < a_cost and b_cost <= mid_cost:
                edge_cost = b_cost
                edge_point = b_point
            elif mid_cost < a_cost and mid_cost < b_cost:
                edge_cost = mid_cost
                edge_point = mid_point

            # If the cost is new
            pathcost = current_dist + edge_cost
            if adj_box not in cost_so_far or pathcost < cost_so_far[adj_box]:
                cost_so_far[adj_box] = pathcost
                came_from[adj_box] = current_box
                box_coord[adj_box] = edge_point
                heappush(priorityQ, (pathcost, adj_box))
    
    if des_box not in came_from.keys():
        return noPathFound(boxes, "Dijkstra's could not reach the destination")

    path = getCoordPathFromDict( came_from, box_coord, source_point, destination_point, src_box, des_box )

    return path, boxes

def a_star( mesh, source_point, destination_point, src_box, des_box ):
    path = []
    boxes = {}

    # back pointer
    came_from = {src_box: None}

    # coordinates of each box pointer
    box_coord = {src_box: source_point}

    # keeps track of cost
    cost_so_far = {src_box: 0}

    priorityQ = [(0, src_box)]

    while priorityQ:
        popped_element = heappop(priorityQ)
        current_box = popped_element[1]
        current_dist = cost_so_far[current_box]

        # Finish the search if we reached the destination
        if current_box == des_box:
            break
        
        # Calculate cost from current box to all the adjacent ones
        for adj_box in mesh['adj'][current_box]:
            boxes[adj_box] = current_box

            # Finds edge range
            x_range = [ max(current_box[0], adj_box[0]), min(current_box[1], adj_box[1]) ]
            y_range = [ max(current_box[2], adj_box[2]), min(current_box[3], adj_box[3]) ]

            # Mid point of edge range
            mid_point = ( (x_range[0] + x_range[1]) * 0.5, (y_range[0] + y_range[1]) * 0.5 )

            # points at the ends of the edge
            a_point = (x_range[0],y_range[0])
            b_point = (x_range[1],y_range[1])

            # Calculates cost and estimate from each edge point
            a_cost = getDistance( box_coord[current_box], a_point )
            a_estimate = a_cost + getDistance( a_point, destination_point )

            b_cost = getDistance( box_coord[current_box], b_point )
            b_estimate = b_cost + getDistance( b_point, destination_point )

            mid_cost = getDistance( box_coord[current_box], mid_point )
            mid_estimate = mid_cost + getDistance( mid_point, destination_point )

            # Lowest cost is added to the queue
            lowest_estimate = min(a_estimate, b_estimate, mid_estimate)
            if lowest_estimate == a_estimate:
                edge_cost = a_cost
                edge_estimate = a_estimate
                edge_point = a_point
            elif lowest_estimate == b_estimate:
                edge_cost = b_cost
                edge_estimate = b_estimate
                edge_point = b_point
            elif lowest_estimate == mid_estimate:
                edge_cost = mid_cost
                edge_estimate = mid_estimate
                edge_point = mid_point

            # If the cost is new
            pathcost = current_dist + edge_cost
            if adj_box not in cost_so_far or pathcost < cost_so_far[adj_box]:
                cost_so_far[adj_box] = pathcost
                came_from[adj_box] = current_box
                box_coord[adj_box] = edge_point
                heappush(priorityQ, (current_dist + edge_estimate, adj_box))
    
    if des_box not in came_from.keys():
        return noPathFound(boxes, "A* could not reach the destination")

    path = getCoordPathFromDict( came_from, box_coord, source_point, destination_point, src_box, des_box )

    return path, boxes

def bidirectional_a_star( mesh, source_point, destination_point, src_box, des_box ):
    path = []
    boxes = {}

    # back pointer
    forward_came_from = {src_box: None}
    backward_came_from = {des_box: None}

    # coordinates of each box pointer
    forward_box_coord = {src_box: source_point}
    backward_box_coord = {des_box: destination_point}

    # keeps track of cost
    forward_cost_so_far = {src_box: 0}
    backward_cost_so_far = {des_box: 0}

    # keys for the goal that is enqueued in the priorityQ
    # to prevent typos
    # if key is destination, then search forward
    # if key is source, then search backwards
    des_key = 'destination'
    src_key = 'source'

    priorityQ = [(0, src_box, des_key), (0, des_box, src_key)]

    path_found = False

    while priorityQ:
        popped_element = heappop(priorityQ)
        current_box = popped_element[1]
        current_goal = popped_element[2]
        if current_goal == des_key:
            current_dist = forward_cost_so_far[current_box]
        else:
            current_dist = backward_cost_so_far[current_box]

        # Finish the search if we reached the destination
        if current_goal == des_key and current_box in backward_came_from:
            path_found = True
            intersection_box = current_box
            break
        elif current_goal == src_key and current_box in forward_came_from:
            path_found = True
            intersection_box = current_box
            break
        
        # Calculate cost from current box to all the adjacent ones
        for adj_box in mesh['adj'][current_box]:
            boxes[adj_box] = current_box

            # Finds edge range
            x_range = [ max(current_box[0], adj_box[0]), min(current_box[1], adj_box[1]) ]
            y_range = [ max(current_box[2], adj_box[2]), min(current_box[3], adj_box[3]) ]

            # Mid point of edge range
            mid_point = ( (x_range[0] + x_range[1]) * 0.5, (y_range[0] + y_range[1]) * 0.5 )

            # points at the ends of the edge
            a_point = (x_range[0],y_range[0])
            b_point = (x_range[1],y_range[1])

            if current_goal == des_key:
                curr_coords = forward_box_coord[current_box]
            else:
                curr_coords = backward_box_coord[current_box]

            # Calculates cost and estimate from each edge point
            a_cost = getDistance( curr_coords, a_point )
            a_estimate = a_cost + getDistance( a_point, destination_point )

            b_cost = getDistance( curr_coords, b_point )
            b_estimate = b_cost + getDistance( b_point, destination_point )

            mid_cost = getDistance( curr_coords, mid_point )
            mid_estimate = mid_cost + getDistance( mid_point, destination_point )

            # Lowest cost is added to the queue
            lowest_estimate = min(a_estimate, b_estimate, mid_estimate)
            if lowest_estimate == a_estimate:
                edge_cost = a_cost
                edge_estimate = a_estimate
                edge_point = a_point
            elif lowest_estimate == b_estimate:
                edge_cost = b_cost
                edge_estimate = b_estimate
                edge_point = b_point
            elif lowest_estimate == mid_estimate:
                edge_cost = mid_cost
                edge_estimate = mid_estimate
                edge_point = mid_point

            # If the cost is new
            pathcost = current_dist + edge_cost
            if current_goal == des_key:
                if adj_box not in forward_cost_so_far or pathcost < forward_cost_so_far[adj_box]:
                    forward_cost_so_far[adj_box] = pathcost
                    forward_came_from[adj_box] = current_box
                    forward_box_coord[adj_box] = edge_point
                    heappush(priorityQ, (current_dist + edge_estimate, adj_box, des_key))
            else:
                if adj_box not in backward_cost_so_far or pathcost < backward_cost_so_far[adj_box]:
                    backward_cost_so_far[adj_box] = pathcost
                    backward_came_from[adj_box] = current_box
                    backward_box_coord[adj_box] = edge_point
                    heappush(priorityQ, (current_dist + edge_estimate, adj_box, src_key))
    
    if not path_found:
        return noPathFound(boxes, "Bidirectional A* could not reach the destination")

    path = getCoordPathForBiDirAStar( forward_came_from, backward_came_from, forward_box_coord, backward_box_coord, source_point, destination_point, intersection_box, src_box, des_box )

    return path, boxes