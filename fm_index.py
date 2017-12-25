from BWT import *


class FMIndex:
    """
        FM Index for string exact matching. Combining BWT with a few auxiliary
        data structure.
    """
    def __init__(self, string, a, b):
        """
        :param string: The string in the database for which to build a FM Index
        :param a: The fraction of rows we keep for the tally
        :param b: The fraction of elements in Suffix Array we keep
        """
        self.string = string
        self.a = a
        self.b = b
        self.bwt = bwtViaBwm(self.string)
        self.fcol = self.generate_column()
        self.tally = self.generate_tally()
        self.sa = self.generate_suffixarray()

    def generate_column(self):
        """
        From the Burrows-Wheeler Transform matrix to get the F column
        :return: F column as a dictionary, the key is the character
        in the string, the value is its number.
        """
        fcol = {}
        bwt_matrix = bwm(self.string)
        for row in bwt_matrix:
            if row[0] in fcol:
                fcol[row[0]] += 1
            elif row[0] != "$":
                fcol[row[0]] = 1
        return fcol

    def generate_tally(self):
        """
        Based on the bwt,the L column, to build a tally table. The
        number of its columns is the number of elements in the character set.
        The number of rows is equal to the length of the string divided by
        the a.
        :return: A tally table as a dictionary. The keys are characters in
        character set, and for each key the value is a list storing the number
        of this character occurring above every ath row.
        """
        tally = {}
        temp = {}
        for e in self.fcol:
            tally[e] = []
            temp[e] = 0
            for i in range(len(self.string) / self.a + 1):
                tally[e].append(0)
        tally[self.bwt[0]][0] = 1
        for i in range(1, len(self.bwt)):
            if self.bwt[i] in temp:
                temp[self.bwt[i]] += 1
            if i % self.a == 0:
                for c in tally:
                    tally[c][i / self.a] = tally[c][i / self.a - 1] + temp[c]
                    temp[c] = 0
        return tally

    def get_rank(self, i, c):
        """
        :param i: The row number of the character c in L column
        :param c: A charcter c
        :return: The B-rank of c in the ith row
        """
        if i % self.a == 0:
            return self.tally[c][i / self.a]
        else:
            count = 0
            while i % self.a != 0:
                if self.bwt[i] == c:
                    count += 1
                i -= 1
            return self.tally[c][i / self.a] + count

    def get_row_num(self, c, r):
        """
        :param c: A character c
        :param r: The B-rank of c
        :return: The row num of the c with rank r in F column.
        """
        result = 0
        for e in self.fcol:
            if e < c:
                result += self.fcol[e]
        result += (r + 1)
        return result

    def bwt_reverse(self):
        """Recover the original string based on F column and L column
        and the tally table
        :return: Original string
        """
        result = ""
        c = self.bwt[0]
        i = 0
        while c != "$":
            rank = self.get_rank(i, c)
            result = c + result
            i = 0
            for e in self.fcol:
                if e < c:
                    i += self.fcol[e]
            i += rank
            c = self.bwt[i]
        return result

    def query(self, qstr):
        """
        To query whether the query string is a substring of the original string.
        :param qstr: Query String
        :return: If the query string is the substring return true. Otherwise return
        false.
        """
        if qstr == "":
            return False, (-1, -1)
        c = qstr[-1]
        if c not in self.fcol:
            return False, (-1, -1)
        srow = 0
        for e in self.fcol:
            if e < c:
                srow += self.fcol[e]
        srow += 1
        erow = srow + self.fcol[c] - 1
        for i in xrange(len(qstr) - 2, -1, -1):
            c = qstr[i]
            if c not in self.fcol:
                return False, (srow, erow)
            srank = self.get_rank(srow - 1, c)
            erank = self.get_rank(erow, c)
            if srank == erank:
                return False, (srow, erow)
            srow = self.get_row_num(c, srank)
            erow = srow + erank - srank - 1
            if erow >= len(self.bwt):
                return False, (-1, -1)
        return True, (srow, erow)

    def generate_suffixarray(self):
        """
        Generate a suffix array with b fraction
        :return: A segment of suffix array
        """
        sa = suffixArray(self.string)
        result = {}
        for i in range(len(self.string)):
            if sa[i] % self.b == 0:
                result[i] = sa[i]
        return result

    def get_offset1(self, qstr):
        """
        Get the offset of the original string using Suffix Array, tally and F column
        :param qstr: Query sequence
        :return: The first matching character's offset in original string
        """
        is_sub, rows = self.query(qstr)
        offsets = []
        if is_sub:
            for i in range(rows[0], rows[1] + 1):
                j = i
                count = 0
                while j not in self.sa:
                    c = self.bwt[j]
                    count += 1
                    rank = self.get_rank(j, c) - 1
                    j = self.get_row_num(c, rank)
                offsets.append(self.sa[j] + count)
        return offsets

    def get_offset(self, qstr, a, b, bwt, fcol, tally, sa):
        """
        Get the offset of the original string using Suffix Array, tally and F column
        :param qstr: Query sequence
        :return: The first matching character's offset in original string
        """

        self.a = a
        self.b = b
        self.bwt = bwt
        self.fcol = fcol
        self.tally = tally
        self.sa = sa


        is_sub, rows = self.query(qstr)
        offsets = []
        if is_sub:
            for i in range(rows[0], rows[1] + 1):
                j = i
                count = 0
                while j not in self.sa:
                    c = self.bwt[j]
                    count += 1
                    rank = self.get_rank(j, c) - 1
                    j = self.get_row_num(c, rank)
                offsets.append(self.sa[j] + count)
        return offsets