def rotations(t):
    """
    :param t: input string
    :return: list of all different rotation of the original string
    """
    tt = t * 2
    return [tt[i:i+len(t)] for i in range(0, len(t))]


def bwm(t):
    """
    :param t: Origianl String
    :return: Burrows-Wheeler Transform matrix
    """
    return sorted(rotations(t))


def bwtViaBwm(t):
    """
    :param t: Original String
    :return:  Last column of the Burrows-Wheeler Transform matrix,
    which is Burrows-Wheeler Transform string
    """
    return ''.join(map(lambda x: x[-1], bwm(t)))


def suffixArray(s):
    """
    :param s: Original string
    :return: A suffix array of the original String
    """
    satups = sorted([(s[i:], i) for i in xrange(0, len(s))])
    return map(lambda x: x[1], satups)





