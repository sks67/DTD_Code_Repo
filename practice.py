#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Learning how to catch errors for the first 
# time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

from __future__ import division

try:
    a = input('a = ')
    b = input('b = ')
    print a/b
except ZeroDivisionError:
    print 'cant divide by 0, bruh'
    print 'a = ', a
    print 'b = ', b
