def fixBadCoordinates(line_parts):
    '''
    Fixes bad coordinates with no spaces between them
    Arguments:
    * line_parts (list): list of bad coordinates
    Returns: (list) fixed x,y,z coordinates for atom
    '''
    coords = []
    chunk = ''.join(line_parts)
    string = ''
    i = 0
    for char in chunk:
        if (char != '.') and (i == 0):
            string = string + char
        if (char != '.') and (i != 0):
            string = string + char
            i += 1
            if i == 4:
                coords.append(float(string))
                if len(coords) == 3:
                    return coords
                i = 0
                string = ''
        if char == '.':
            string = string + char
            i += 1