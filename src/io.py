"""Module for performing input/output operations.

Functions are:

*ParseInfile* : reads key / value pairs from an input file.

*ParseBoolValue* : reads Boolean values from input dictionary.

*ParseIntValue* : reads integer values from input dictionary.

*ParseFloatValue* : reads float values from input dictionary.

*ParseStringValue* : reads string values from input dictionary.

*ParseSeqValue* : reads sequence from file specified by input dictionary.

*ParseFileList* : reads file list specified by input dictionary.


Written by Jesse Bloom.
"""


import os
import sys
import re
import cStringIO


def ParseInfile(f):
    """Reads key / value pairs from an input file.

    *f* should be a readable file-like object.

    Starting at the current position in *f*, reads all remaining lines
    until the end of the file-like object. Any blank line or line
    with a first character of # is disregarded (# indicates a comment
    line). All other lines should contain exactly two entries, the first
    being the key and the second being the value. The key is construed
    to be the first word, and cannot have any spaces. The value is all
    of the text that follows the key up to the newline. 
    The key / value pairs are returned in a
    dictionary, as illustrated in the example below. If there is a
    duplicate key, and exception is raised.

    Example of successfully reading two key/value pairs

    >>> f = cStringIO.StringIO()
    >>> f.write('# comment line followed by blank line\\n\\n')
    >>> f.write('key1 first_value\\n')
    >>> f.write('# now another key with two values\\nkey2 2 value2')
    >>> f.seek(0)
    >>> ParseInfile(f) == {'key1':'first_value', 'key2':'2 value2'}
    True

    Example of duplicate key name

    >>> f = cStringIO.StringIO()
    >>> f.write('# comment line followed by blank line\\n\\n')
    >>> f.write('key1 first_value\\n')
    >>> f.write('# now another key\\nkey2 2')
    >>> f.write('\\nkey1 1')
    >>> f.seek(0)
    >>> ParseInfile(f) == {'key1':'first_value', 'key2':'2'}
    Traceback (most recent call last):
        ...
    ValueError: duplicate key: key1

    """
    d = {}
    for line in f:
        line = line.strip()
        if line and (not (line.isspace() or line[0] == '#')):
            entries = line.split(None, 1)
            if len(entries) != 2:
                raise ValueError("invalid value/key, not two entries:"\
                        + ' %s' % line.strip())
            key = entries[0].strip()
            value = entries[1].strip()
            if key in d:
                raise ValueError("duplicate key: %s" % key)
            d[key] = value
    return d


def ParseFileList(d, key):
    """Gets list of files from input dictionary.

    *d* is a dictionary of key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that
    should give a list of one or more file names.

    Returns a list of the file names specified by *key*. If one or more of
    these files does not exist, raises and IOError. If key does not
    exist in *d*, raises and KeyError.
    """
    if key not in d:
        raise KeyError("Did not find key of %s" % key)
    files = []
    for f in d[key].split():
        f = f.strip()
        if not os.path.isfile(f):
            raise IOError("Cannot find specified file %s" % f)
        files.append(f)
    return files


def ParseBoolValue(d, key):
    """Gets Boolean argument from input dictionary.

    *d* is a dictionary of key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that
    should be *True* or *False*.

    Returns the Boolean truth value specified by *key*. Raises a *ValueError*
    if *key* does not exist in *d*, and a *ValueError* if it specifies
    something other than *True* or *False*.

    >>> d = {'gzipped':'True', 'applyfilter':'False', 'a1':'a1.txt'}
    >>> ParseBoolValue(d, 'gzipped')
    True

    >>> ParseBoolValue(d, 'applyfilter')
    False

    >>> ParseBoolValue(d, 'a1')
    Traceback (most recent call last):
       ...
    ValueError: value a1.txt for a1 is not True/False

    >>> ParseBoolValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    elif d[key] == 'True':
        return True
    elif d[key] == 'False':
        return False
    else:
        raise ValueError("value %s for %s is not True/False" % (d[key],\
                key))


def ParseIntValue(d, key):
    """Gets integer argument from input dictionary.

    *d* is a dictionary of key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that
    should be an integer.

    Returns the integer specified by *key*. Raises a *ValueError*
    if *key* does not exist in *d*, and a *ValueError* if it specifies
    something other than an integer.

    >>> d = {'maxn':'2', 'minq':'20.5'}
    >>> ParseIntValue(d, 'maxn')
    2

    >>> ParseIntValue(d, 'minq')
    Traceback (most recent call last):
       ...
    ValueError: value 20.5 for minq is not integer

    >>> ParseIntValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    try:
        return int(d[key])
    except ValueError:
        raise ValueError('value %s for %s is not integer' % \
                (d[key], key))


def ParseFloatValue(d, key):
    """Gets float argument from input dictionary.

    *d* a dictionary of key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that
    should be convertible to a float.

    Returns the float specified by *key*. Raises a *ValueError*
    if *key* does not exist in *d*, and a *ValueError* if it specifies
    something other than an integer.

    >>> d = {'maxn':'2', 'minq':'20.5', 'string':'hello'}
    >>> print "%.1f" % ParseFloatValue(d, 'maxn')
    2.0

    >>> print "%.1f" % ParseFloatValue(d, 'minq')
    20.5

    >>> ParseFloatValue(d, 'string')
    Traceback (most recent call last):
       ...
    ValueError: value hello for string is not float

    >>> ParseFloatValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError('did not find a key for %s' % key)
    try:
        return float(d[key])
    except ValueError:
        raise ValueError('value %s for %s is not float' % \
                (d[key], key))


def ParseSeqValue(d, key):
    """Reads sequence from FASTA file specified by input dictionary.

    *d* is a dictionary of key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that is a
    filename of a readable file that contains exactly one sequence
    in FASTA format.

    Returns a string corresponding to the sequence specified in the
    FASTA file. Raises a *ValueError* if the filename does not exist, or
    does not specify exactly one sequence.
    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    filename = d[key]
    if not os.path.isfile(filename):
        raise ValueError("Cannot find file %s specified by %s" % \
                (filename, key))
    seqs = mapmuts.sequtils.ReadFASTA(filename)
    if len(seqs) != 1:
        raise ValueError("file %s does not contain exactly one"\
                % filename + ' sequence')
    return seqs[0][1]


def ParseStringValue(d, key):
    """Reads string argument specified by input dictionary.

    *d* is a dictionary of the key/value string pairs as returned by
    **ParseInfile**.

    *key* is a string in this dictionary that specifies a value that
    is a string.

    Returns string corresponding to *key*. Raises a *ValueError* if 
    the *key* does not exist in *d*.

    >>> d = {'outfileprefix':'test_example', 'samplename':'test example'}
    >>> ParseStringValue(d, 'outfileprefix')
    'test_example'

    >>> ParseStringValue(d, 'samplename')
    'test example'

    >>> ParseStringValue(d, 'otherkey')
    Traceback (most recent call last):
       ...
    ValueError: did not find a key for otherkey

    """
    if key not in d:
        raise ValueError("did not find a key for %s" % key)
    return d[key].strip()




if __name__ == '__main__':
    import doctest
    doctest.testmod()
