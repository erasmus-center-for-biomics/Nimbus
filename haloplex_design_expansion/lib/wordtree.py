#!/bin/env python

__author__ = 'rwwbr_000'

import sys

'''

This file contains the WordTree class. The CMPTree holds words in a tree like structure
and (should) allows for rapid querying. Some pre processing is involved.
'''


class CharNode(object):

    def __init__(self, c='^' ):
        self._char = c
        self._children = []
        self._parent = None

    def add_child(self, oth):
        for c in self._children:
            if c == oth:
                return c
        oth._parent = self
        self._children.append(oth)
        return oth

    def __eq__(self, other):
        if self._char == other._char :
            return True
        else:
            return False

    def parent(self):
        self._parent

    def is_root(self):
        if self._parent == None:
            return True
        else:
            return False

class WordTree(object):

    def __init__(self):
        self._root        = CharNode( c=None )
        self._words_added = 0

        # a special character that will reset the protected
        # zone in the matching algorithm
        self._reset_char  = '|'
        self._protected   = 5
        self._maxscore    = 1

    def add_word(self, word=""):
        '''
        Adds a new word to the tree
        :param word: the word to add
        :return: nothing
        '''
        cur = self._root
        word = word + '$'
        for i in xrange(len(word)):
            c   = word[i]
            x   = CharNode( c=c )
            cur = cur.add_child( x )

        self._words_added += 1

    def is_present( self, word="", safezone=None, node=None, score=0 ):
        '''
        Determine whether the sequence is present in the tree

        conditions on when the loop terminates are:
            1) word length is finished
            - per child:
                2) end node in tree is reached
                3) score exceeds maximum score
                4) downstream hit is found

        :param word: the word to match
        :param node: the current node
        :param score: the score thus far
        :return: present True or False
        '''
        # the default value of the function is False
        rval = False

        # starting condition
        if node is None:
            node = self._root

        #
        if safezone is None:
            safezone = self._protected

        # we reached the end of the line and are done
        if len(word) == 0:
            rval = True
        else:

            # reset the scores if we find the reset character
            if word[0] == self._reset_char:
                safezone = self._protected + 1 # add 1 to move over this character
                score    = 0

            # foreach child in the node
            for child in node._children:

                # set the temp score for the current node
                tscore = score
                # sys.stderr.write( "%d\t%d\t%s\t%s\t%d\n" %(tscore, safezone, word, child._char, len(node._children)) )

                # are we at the end of the index?
                if child._char == '$':
                    rval = True
                    break
                # if we are still dealing with protected bases
                # make sure we don't have a mismatch, otherwise next iteration
                elif safezone >= 0 and child._char != word[0] :
                    continue
                # or should we increase the score
                elif child._char != word[0]:
                    tscore += 1

                # terminate when the score exceed our maximum
                if tscore <= self._maxscore:

                    # enter the next cycle
                    rval = self.is_present( word[1:], node=child, safezone=safezone-1, score=tscore )
                    if rval == True:
                        break
                else:
                    continue

        # when all else fails
        return rval
