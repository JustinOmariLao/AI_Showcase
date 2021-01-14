from p1_support import load_level, show_level, save_level_costs
from math import inf, sqrt
from heapq import heappop, heappush, heapify
import sys


def dijkstras_shortest_path(initial_position, destination, graph, adj):
    """ Searches for a minimal cost path through a graph using Dijkstra's algorithm.

    Args:
        initial_position: The initial cell from which the path extends.
        destination: The end location for the path.
        graph: A loaded level, containing walls, spaces, and waypoints.
        adj: An adjacency function returning cells adjacent to a given cell as well as their respective edge costs.

    Returns:
        If a path exits, return a list containing all cells from initial_position to destination.
        Otherwise, return None.

    """

    # assembling queue on heap
    priorityQ = []
    heappush(priorityQ, (0, initial_position))

    # back pointer
    came_from = {initial_position: None}

    # keeps track of cost
    cost_so_far = {initial_position: 0}

    while priorityQ:
        current_cost, current_node = heappop(priorityQ)

        if current_node == destination:
            print("Total cost = ", cost_so_far[current_node], '\n')
            curr = current_node
            path = []
            while curr != None:
                path.append(curr)
                curr = came_from[curr]
            return path

        for new_node, new_cost in navigation_edges(graph, current_node):
            pathCost = new_cost + current_cost
            if new_node not in cost_so_far or pathCost < cost_so_far[new_node]:
                cost_so_far[new_node] = pathCost
                heappush(priorityQ, (pathCost, new_node))
                came_from[new_node] = current_node
    return None


def find_node(node, pqueue):
    for element in pqueue:
        if element[1] == node:
            return True
    return False


def dijkstras_shortest_path_to_all(initial_position, graph, adj):
    """ Calculates the minimum cost to every reachable cell in a graph from the initial_position.

    Args:
        initial_position: The initial cell from which the path extends.
        graph: A loaded level, containing walls, spaces, and waypoints.
        adj: An adjacency function returning cells adjacent to a given cell as well as their respective edge costs.

    Returns:
        A dictionary, mapping destination cells to the cost of a path from the initial_position.
    """

    # assembling queue on heap
    priorityQ = []
    heappush(priorityQ, (0, initial_position))

    # back pointer
    came_from = {initial_position: None}

    # keeps track of cost
    cost_so_far = {initial_position: 0}

    while priorityQ:
        current_cost, current_node = heappop(priorityQ)

        for new_node, new_cost in navigation_edges(graph, current_node):
            pathCost = new_cost + current_cost
            if new_node not in cost_so_far or pathCost < cost_so_far[new_node]:
                cost_so_far[new_node] = pathCost
                heappush(priorityQ, (pathCost, new_node))
                came_from[new_node] = current_node

    return cost_so_far

    pass


def navigation_edges(level, cell):
    """ Provides a list of adjacent cells and their respective costs from the given cell.

    Args:
        level: A loaded level, containing walls, spaces, and waypoints.
        cell: A target location.

    Returns:
        A list of tuples containing an adjacent cell's coordinates and the cost of the edge joining it and the
        originating cell.

        E.g. from (0,0):
            [((0,1), 1),
             ((1,0), 1),
             ((1,1), 1.4142135623730951),
             ... ]
    """
    if cell in level['walls']:  # skip this call
        return []

    x = cell[0]  # x = cell x-coordinate
    y = cell[1]  # x = cell y-coordinate

    # list that holds adjacencies and costs
    adj = []

    # iterates in around cell
    for a in range(x - 1, x + 2, 1):  # traverses the rows of the cell
        for b in range(y - 1, y + 2, 1):  # traverses the columns of the cell
            if (a, b) in level['walls']:  # if wall then skip iteration
                continue
            elif (a, b) == cell:  # if original cell then skip iteration
                continue
            else:
                if (a, b) == (x - 1, y + 1) or (a, b) == (x + 1, y + 1) or (a, b) == (x - 1, y - 1) or (a, b) == (
                        x + 1, y - 1):  # if diagonal then use distance formula
                    dCost = (0.5 * sqrt(2) * level['spaces'][(a, b)]) + (0.5 * sqrt(2) * level['spaces'][(x, y)])
                    adj.append(((a, b), dCost))
                else:  # if horizontal/vertical space then continue
                    sCost = (0.5 * level['spaces'][(a, b)]) + (0.5 * level['spaces'][(x, y)])
                    adj.append(((a, b), sCost))

    return adj


def test_route(filename, src_waypoint, dst_waypoint):
    """ Loads a level, searches for a path between the given waypoints, and displays the result.

    Args:
        filename: The name of the text file containing the level.
        src_waypoint: The character associated with the initial waypoint.
        dst_waypoint: The character associated with the destination waypoint.

    """

    # Load and display the level.
    level = load_level(filename)
    show_level(level)

    # Retrieve the source and destination coordinates from the level.
    src = level['waypoints'][src_waypoint]
    dst = level['waypoints'][dst_waypoint]

    # Search for and display the path from src to dst.
    path = dijkstras_shortest_path(src, dst, level, navigation_edges)
    if path:
        show_level(level, path)
    else:
        print("No path possible!")


def cost_to_all_cells(filename, src_waypoint, output_filename):
    """ Loads a level, calculates the cost to all reachable cells from
    src_waypoint, then saves the result in a csv file with name output_filename.

    Args:
        filename: The name of the text file containing the level.
        src_waypoint: The character associated with the initial waypoint.
        output_filename: The filename for the output csv file.

    """

    # Load and display the level.
    level = load_level(filename)
    show_level(level)

    # Retrieve the source coordinates from the level.
    src = level['waypoints'][src_waypoint]

    # Calculate the cost to all reachable cells from src and save to a csv file.
    costs_to_all_cells = dijkstras_shortest_path_to_all(src, level, navigation_edges)
    save_level_costs(level, costs_to_all_cells, output_filename)



if __name__ == '__main__':
    sys.stdout = open('test_maze_path.txt', 'w')
    test_route('test_maze.txt', 'a', 'd')

    sys.stdout = open('console_output.txt', 'w')

    test_route('my_maze.txt', 'a', 'b')
    test_route('my_maze.txt', 'a', 'c')
    test_route('my_maze.txt', 'a', 'd')
    test_route('my_maze.txt', 'a', 'e')

    # Use this function to calculate the cost to all reachable cells from an origin point.
    cost_to_all_cells('my_maze.txt', 'a', 'my_maze_costs.csv')
