# Made by
#   Justin Lao, jlao3
#   Claudio Sangeroki, csangero

from mcts_node import MCTSNode
from random import choice
from math import sqrt, log, inf
import p3_t3

num_nodes = 1000
num_nodesTWO = 1000
explore_faction = 2.


def traverse_nodes(node, board, state, identity):
    """ Traverses the tree until the end criterion are met.

    Args:
        node:       A tree node from which the search is traversing.
        board:      The game setup.
        state:      The state of the game.
        identity:   The bot's identity, either 'red' or 'blue'.

    Returns:        A node from which the next stage of the search can proceed.

    """
    currentNode = node

    if currentNode.untried_actions or not currentNode.child_nodes:
        return currentNode
    else:
        nextNode = None
        highScore = -inf
        # catcher contains values we don't care about in node, childNode is what we want to check
        for catcher, childNode in currentNode.child_nodes.items():
            childNode.visits += 1
            if identity == 1:   # if red then value is equal to 1, thus nothing happens
                tempScore = (childNode.wins / childNode.visits) + (explore_faction *
                                                                   sqrt(log(childNode.parent.visits) / childNode.visits))
            else:   # else identity is blue, thus multiple by -1
                tempScore = ((-1 * childNode.wins) / childNode.visits) + (
                    explore_faction * sqrt(log(childNode.parent.visits) / childNode.visits))
            if tempScore > highScore:   # if temporary score is higher then our best, then make that the next searched node and update highScore
                highScore = tempScore
                nextNode = childNode
        return traverse_nodes(nextNode, board, state, identity)
    # Hint: return leaf_node


def expand_leaf(node, board, state):
    """ Adds a new leaf to the tree by creating a new child node for the given node.

    Args:
        node:   The node for which a child will be added.
        board:  The game setup.
        state:  The state of the game.

    Returns:    The added child node.

    """
    try:
        action = node.untried_actions.pop()  # tries to pop() an action from the list
    except:
        return node  # if list is empty, return node
    else:
        state = board.next_state(state, action)
        new_node = MCTSNode(parent=node, parent_action=action, action_list=board.legal_actions(
            state))  # create new leaf with list of legal actions
        # make pointer of child node equal the new_node
        node.child_nodes[action] = new_node
        return new_node


def rollout(board, state):
    """ Given the state of the game, the rollout plays out the remainder randomly.

    Args:
        board:  The game setup.
        state:  The state of the game.

    """
    while not board.is_ended(state):    # while game is not done
        possible_actions = board.legal_actions(
            state)   # acquire a list of legal actions
        action = choice(possible_actions)   # choose a random action
        # state becomes the next state of chosen action
        state = board.next_state(state, action)
    return state


def backpropagate(node, won):
    """ Navigates the tree from a leaf node to the root, updating the win and visit count of each node along the path.

    Args:
        node:   A leaf node.
        won:    An indicator of whether the bot won or lost the game.

    """
    if node != None:    # while backtracking is true
        node.wins += won    # -1 if lost, 0 for draw/nothing, 1 for win
        node.visits += 1    # count for visits
        node = node.parent  # traverse to parent node
        backpropagate(node, won)


def think(board, state):
    """ Performs MCTS by sampling games and calling the appropriate functions to construct the game tree.

    Args:
        board:  The game setup.
        state:  The state of the game.

    Returns:    The action to be taken.

    """
    identity_of_bot = board.current_player(state)
    root_node = MCTSNode(parent=None, parent_action=None,
                         action_list=board.legal_actions(state))

    if identity_of_bot == 1:
        for steps in range(num_nodes):
            # Copy the game for sampling a playthrough
            sampled_game = state

            # Start at root
            node = root_node

            # Do MCTS - This is all you!

            # traverse tree until we find best leaf and select it
            node = traverse_nodes(node, board, sampled_game, identity_of_bot)

            # creating a list of actions done by leaf node
            select_actions = aList(node)

            # apply the next state for each action in our list
            for action in select_actions:
                sampled_game = board.next_state(sampled_game, action)

            # if there are no more untried actions, determine winner by point value
            if not node.untried_actions:
                won = board.points_values(sampled_game)[1]
            else:   # else add leaf node, sample next state with the parent of new node and rollout
                node = expand_leaf(node, board, sampled_game)
                sampled_game = board.next_state(
                    sampled_game, node.parent_action)
                sampled_game = rollout(board, sampled_game)
                won = board.points_values(sampled_game)[1]  # determine winner
            backpropagate(node, won)    # update wins and visits of node

        # determine best action, works the same as traverse_nodes
        best_action = bestAction(root_node, identity_of_bot)
        return best_action
    else:
        for steps in range(num_nodesTWO):
            # Copy the game for sampling a playthrough
            sampled_game = state

            # Start at root
            node = root_node

            # Do MCTS - This is all you!

            # traverse tree until we find best leaf and select it
            node = traverse_nodes(node, board, sampled_game, identity_of_bot)

            # creating a list of actions done by leaf node
            select_actions = aList(node)

            # apply the next state for each action in our list
            for action in select_actions:
                sampled_game = board.next_state(sampled_game, action)

            # if there are no more untried actions, determine winner by point value
            if not node.untried_actions:
                won = board.points_values(sampled_game)[1]
            else:   # else add leaf node, sample next state with the parent of new node and rollout
                node = expand_leaf(node, board, sampled_game)
                sampled_game = board.next_state(
                    sampled_game, node.parent_action)
                sampled_game = rollout(board, sampled_game)
                won = board.points_values(sampled_game)[1]  # determine winner
            backpropagate(node, won)    # update wins and visits of node

        # determine best action, works the same as traverse_nodes
        best_action = bestAction(root_node, identity_of_bot)
        return best_action


def aList(node):
    # creating a list of actions done by leaf node
    selected_node = node
    select_actions = []
    while selected_node.parent:
        select_actions.append(selected_node.parent_action)
        selected_node = selected_node.parent
    select_actions.reverse()
    return select_actions


def bestAction(root_node, identity_of_bot):
    best_winrate = -inf
    if identity_of_bot == 1:
        sign = 1
    else:
        sign = -1
    for action, child in root_node.child_nodes.items():
        child_winrate = (child.wins / child.visits) * sign
        if child_winrate > best_winrate:
            best_action = action
            best_winrate = child_winrate

    return best_action
