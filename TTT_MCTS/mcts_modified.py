# Made by
#   Justin Lao, jlao3
#   Claudio Sangeroki, csangero

from mcts_node import MCTSNode
from random import choice
from math import sqrt, log, inf
import random
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
    # Winning States

    p = board.current_player(state)

    while not board.is_ended(state):    # while game is not done
        winConOne = False
        winConTwo = False
        # Evaluate Score 1000f1 + 100f2 - 10f3 - 1f2
        possible_actions = board.legal_actions(state)
        possible_actions2 = board.legal_actions(state)
        possible_actions3 = board.legal_actions(state)
        movesDone = board.owned_boxes(state)

        sampleMoves = random.choices(possible_actions, k=3)

        for i in sampleMoves:
            movesDone[tuple([i[2], i[3]])
                      ] = p
            if ((movesDone[(0, 0)] == p and movesDone[(0, 1)] == p and movesDone[(0, 2)] == p) or
                (movesDone[(1, 0)] == p and movesDone[(1, 1)] == p and movesDone[(1, 2)] == p) or
                (movesDone[(2, 0)] == p and movesDone[(2, 1)] == p and movesDone[(2, 2)] == p) or
                (movesDone[(0, 0)] == p and movesDone[(1, 0)] == p and movesDone[(2, 0)] == p) or
                (movesDone[(0, 1)] == p and movesDone[(1, 1)] == p and movesDone[(2, 1)] == p) or
                (movesDone[(0, 2)] == p and movesDone[(1, 2)] == p and movesDone[(2, 2)] == p) or
                (movesDone[(0, 0)] == p and movesDone[(1, 1)] == p and movesDone[(2, 2)] == p) or
                    (movesDone[(2, 0)] == p and movesDone[(1, 1)] == p and movesDone[(0, 2)] == p)):
                winConOne = True
                action = i
            movesDone = board.owned_boxes(state)

        if winConOne == True:
            state = board.next_state(state, action)

        # while winConTwo == False and possible_actions2:
        for i in sampleMoves:

            movesDone[tuple([i[2], i[3]])
                      ] = p
            if (((movesDone[(0, 0)] == p) and (movesDone[(0, 1)] != p and movesDone[(0, 1)] != 0) and (movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(0, 1)] == p) and (movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(0, 1)] != p and movesDone[(0, 1)] != 0) and (movesDone[(0, 2)] == p)) or
                ((movesDone[(1, 0)] == p) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(1, 2)] != p and movesDone[(1, 2)] != 0)) or
                ((movesDone[(1, 0)] != p and movesDone[(1, 0)] != 0) and (movesDone[(1, 1)] == p) and (movesDone[(1, 2)] != p and movesDone[(1, 2)] != 0)) or
                ((movesDone[(1, 0)] != p and movesDone[(1, 0)] != 0) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(1, 2)] == p)) or
                ((movesDone[(2, 0)] == p) and (movesDone[(2, 1)] != p and movesDone[(2, 1)] != 0) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0) and (movesDone[(2, 1)] == p) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0) and (movesDone[(2, 1)] != p and movesDone[(2, 1)] != 0) and (movesDone[(2, 2)] == p)) or
                ((movesDone[(0, 0)] == p) and (movesDone[(1, 0)] != p and movesDone[(1, 0)] != 0) and (movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(1, 0)] == p) and (movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(1, 0)] != p and movesDone[(1, 0)] != 0) and (movesDone[(2, 0)] == p)) or
                ((movesDone[(0, 1)] == p) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 1)] != p and movesDone[(2, 1)] != 0)) or
                ((movesDone[(0, 1)] != p and movesDone[(0, 1)] != 0) and (movesDone[(1, 1)] == p) and (movesDone[(2, 1)] != p and movesDone[(2, 1)] != 0)) or
                ((movesDone[(0, 1)] != p and movesDone[(0, 1)] != 0) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 1)] == p)) or
                ((movesDone[(0, 2)] == p) and (movesDone[(1, 2)] != p and movesDone[(1, 2)] != 0) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0) and (movesDone[(1, 2)] == p) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0) and (movesDone[(1, 2)] != p and movesDone[(1, 2)] != 0) and (movesDone[(2, 2)] == p)) or
                ((movesDone[(0, 0)] == p) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(1, 1)] == p) and (movesDone[(2, 2)] != p and movesDone[(2, 2)] != 0)) or
                ((movesDone[(0, 0)] != p and movesDone[(0, 0)] != 0) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 2)] == p)) or
                ((movesDone[(0, 2)] == p) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0)) or
                ((movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0) and (movesDone[(1, 1)] == p) and (movesDone[(2, 0)] != p and movesDone[(2, 0)] != 0)) or
                    ((movesDone[(0, 2)] != p and movesDone[(0, 2)] != 0) and (movesDone[(1, 1)] != p and movesDone[(1, 1)] != 0) and (movesDone[(2, 0)] == p))):
                winConTwo = True
                action = i
            movesDone = board.owned_boxes(state)

        # count = 0
        # if p == 1:
        #     while winConTwo == False and possible_actions2 and count < 10:
        #         action = possible_actions2.pop()
        #         movesDone[tuple([action[2], action[3]])
        #                   ] = board.current_player(state)
        #         if ((movesDone[(0, 0)] == 1 and movesDone[(0, 1)] == 2 and movesDone[(0, 2)] == 2) or
        #             (movesDone[(0, 0)] == 2 and movesDone[(0, 1)] == 1 and movesDone[(0, 2)] == 2) or
        #             (movesDone[(0, 0)] == 2 and movesDone[(0, 1)] == 2 and movesDone[(0, 2)] == 1) or
        #                 (movesDone[(1, 0)] == 1 and movesDone[(1, 1)] == 2 and movesDone[(1, 2)] == 2) or
        #                 (movesDone[(1, 0)] == 2 and movesDone[(1, 1)] == 1 and movesDone[(1, 2)] == 2) or
        #                 (movesDone[(1, 0)] == 2 and movesDone[(1, 1)] == 2 and movesDone[(1, 2)] == 1) or
        #                 (movesDone[(2, 0)] == 1 and movesDone[(2, 1)] == 2 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(2, 0)] == 2 and movesDone[(2, 1)] == 1 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(2, 0)] == 2 and movesDone[(2, 1)] == 2 and movesDone[(2, 2)] == 1) or
        #                 (movesDone[(0, 0)] == 1 and movesDone[(1, 0)] == 2 and movesDone[(2, 0)] == 2) or
        #                 (movesDone[(0, 0)] == 2 and movesDone[(1, 0)] == 1 and movesDone[(2, 0)] == 2) or
        #                 (movesDone[(0, 0)] == 2 and movesDone[(1, 0)] == 2 and movesDone[(2, 0)] == 1) or
        #                 (movesDone[(0, 1)] == 1 and movesDone[(1, 1)] == 2 and movesDone[(2, 1)] == 2) or
        #                 (movesDone[(0, 1)] == 2 and movesDone[(1, 1)] == 1 and movesDone[(2, 1)] == 2) or
        #                 (movesDone[(0, 1)] == 2 and movesDone[(1, 1)] == 2 and movesDone[(2, 1)] == 1) or
        #                 (movesDone[(0, 2)] == 1 and movesDone[(1, 2)] == 2 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(0, 2)] == 2 and movesDone[(1, 2)] == 1 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(0, 2)] == 2 and movesDone[(1, 2)] == 2 and movesDone[(2, 2)] == 1) or
        #                 (movesDone[(0, 0)] == 1 and movesDone[(1, 1)] == 2 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(0, 0)] == 2 and movesDone[(1, 1)] == 1 and movesDone[(2, 2)] == 2) or
        #                 (movesDone[(0, 0)] == 2 and movesDone[(1, 1)] == 2 and movesDone[(2, 2)] == 1) or
        #                 (movesDone[(2, 0)] == 1 and movesDone[(1, 1)] == 2 and movesDone[(0, 2)] == 2) or
        #                 (movesDone[(2, 0)] == 2 and movesDone[(1, 1)] == 1 and movesDone[(0, 2)] == 2) or
        #                 (movesDone[(2, 0)] == 2 and movesDone[(1, 1)] == 2 and movesDone[(0, 2)] == 1)):
        #             winConTwo = True

        #         movesDone = board.owned_boxes(state)
        #         count += 1

        if winConTwo == True and winConOne != True:
            state = board.next_state(state, action)

        else:
            action = choice(possible_actions3)
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

# def evaluate_score(board, state):

#     identity_of_bot = board.current_player(state)

#     localBoardState = board.owned_boxes(state)

#     possible_actions = board.legal_actions(state)

#     scoreDict = {(0, 0): 3, (0, 1): 2, (0, 2): 3, (1, 0): 2, (1, 1)
#                   : 4, (1, 2): 2, (2, 0): 3, (2, 1): 2, (2, 2): 3}
#     # Take middle value if board is empty

#     neighbors = {}
#     neighborList = []

#     for move in possible_actions:
#         for i in range(move[2]-1, move[2]+2):
#             for j in range(move[3]-1, move[3]+2):
#                 if (i != move[2] and j != move[3]):
#                     validMove = tuple([i, j])
#                     if (validMove in localBoardState) and localBoardState[validMove] == 0:
#                         neighborList.append(
#                             tuple([validMove[0], validMove[1]]))
#                 elif (i != move[2] or j != move[3]):
#                     validMove = tuple([i, j])
#                     if (validMove in localBoardState) and localBoardState[validMove] == 0:
#                         neighborList.append(
#                             tuple([validMove[0], validMove[1]]))

#         neighbors[move] = neighborList
#         neighborList = []

#     for move in possible_actions:

#         for i in range(len(neighbors[move])):
#             if localBoardState[neighbors[move][i]] == identity_of_bot:
#                 scoreDict[tuple([move[2], move[3]])] += 10
#             elif localBoardState[neighbors[move][i]] != identity_of_bot and not 0:
#                 scoreDict[tuple([move[2], move[3]])] -= 10
#             else:
#                 scoreDict[tuple([move[2], move[3]])] += 0

#     maxKey = max(scoreDict, key=scoreDict.get)
#     tempMove2 = possible_actions[0]
#     maxMove = tuple([tempMove2[0], tempMove2[1], maxKey[0], maxKey[1]])
#     return maxMove
