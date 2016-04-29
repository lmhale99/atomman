from DataModelDict import DataModelDict
import numpy as np
import atomman as am
import atomman.unitconvert as uc
import sys

def load(cif_file, data_set = None):
    """Read in a cif file or DataModelDict of a cif file and return a System."""    
    try:
        alldict = model(cif_file)
    except:
        alldict = DataModelDict(cif_file)
    
    if data_set is None:
        data_set = -1
    if isinstance(data_set, (int, long)):
        data_set = alldict.keys()[data_set]
    if len(data_set) < 5 or data_set[:5] != 'data_':
        data_set = 'data_' + data_set
    
    alldict = alldict[data_set]
    
    #Read box cell parameters
    try:
        a =     __tofloat(alldict.get('_cell_length_a') )
        b =     __tofloat(alldict.get('_cell_length_b') )
        c =     __tofloat(alldict.get('_cell_length_c') )
        alpha = __tofloat(alldict.get('_cell_angle_alpha') )
        beta =  __tofloat(alldict.get('_cell_angle_beta') )
        gamma = __tofloat(alldict.get('_cell_angle_gamma') )
    except:
        raise ValueError('Invalid cell parameters')
    
    a = uc.set_in_units(a, 'angstrom')
    b = uc.set_in_units(b, 'angstrom')
    c = uc.set_in_units(c, 'angstrom')        
    box = am.Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
          
    #Read atom site fractions
    xlist = alldict.aslist('_atom_site_fract_x')
    ylist = alldict.aslist('_atom_site_fract_y')
    zlist = alldict.aslist('_atom_site_fract_z')
 
    if len(xlist) == 0 or len(xlist) != len(ylist) or len(xlist) != len(zlist):
        raise ValueError('Invalid atom site fractions')
    try:
        site_fracts = []
        for i in xrange(len(xlist)):
            site_fracts.append(np.array([__tofloat(xlist[i]),
                                         __tofloat(ylist[i]), 
                                         __tofloat(zlist[i])] ))    
    except:
        raise ValueError('Invalid atom site fractions')
        
    #Read in symmetry lists
    symms = alldict.aslist('_space_group_symop_operation_xyz')
    if len(symms) == 0:
        symms = alldict.aslist('_symmetry_equiv_pos_as_xyz')
    if len(symms) == 0:
        raise ValueError('No symmetries listed')
    
    for i in xrange(len(symms)):
        if len(symms[i]) > 2:
            if (symms[i][0] == '"' and symms[i][-1] == '"') or (symms[i][0] == "'" and symms[i][-1] == "'"):
                symms[i] = symms[i][1:-1] 
        symms[i] = symms[i].lower()
        symms[i] = symms[i].split(',')
        if len(symms[i]) != 3:
            raise ValueError('Bad symmetry terms')

    site = 0
    atype = []
    pos = []
    for site_fract in site_fracts:
        coords = []
        site += 1

        for symm in symms:
            coord = __calc(symm, site_fract)
            unique = True
            for c in coords:
                if np.allclose(c, coord, rtol = 0.000001):
                    unique = False
                    break
            if unique:
                coords.append(coord)
        for coord in coords:
            atype.append(site)
            pos.append(coord)
    atoms = am.Atoms(natoms=len(atype), prop={'atype':atype, 'pos':pos})
    
    system = am.System(box=box, atoms=atoms, scale=True)
       
    elements = alldict.aslist('_atom_site_type_symbol')
    if len(elements) == 0:
        elements = alldict.aslist('_atom_site_label')
    if len(elements) == 0:
        raise ValueError('No atom site symbols or labels!')
        
    for i in xrange(len(elements)):
        for j in xrange(len(elements[i])):
            if elements[i][j].isdigit() or elements[i][j] in '(+-':
                elements[i] = elements[i][:j]
                break
        if len(elements[i]) == 0:
            raise ValueError('Bad atom site symbol or label!')
    
    return system, elements

def model(cif_file):
    """Read in a cif file and return a DataModelDict of the data."""
    if hasattr(cif_file, 'read'):
        cif_file = cif_file.read()
    alldict = __parse(cif_file)
    
    return alldict
    
def __parse(model):
    """Parses a cif file to an equivalent DataModelDict."""
    lines = model.split('\n')     
    terms = []
    
    #Strip away comments
    i = 0
    while i < len(lines):
        if len(lines[i]) == 0:
            i += 1
        
        #strip comments
        elif '#' in lines[i][0]:
            hash_index = lines[i].index('#')
            lines[i] = lines[i][:hash_index]

        
        #Handle semicolon deliminters 
        elif lines[i][0] == ';':
            term = lines[i][1:]
            i += 1
            while lines[i].strip() != ';':
                term += lines[i]
                i += 1
                if i == len(lines):
                    raise ValueError('Bad file!')
            terms.append(term.strip())
            i += 1
            
        else:
            words = lines[i].split()
    
            #Separate out terms using whitespace
            j = 0
            while j < len(words):
                #Handle double quotes
                if words[j][0] == '"':
                    words[j] = words[j][1:]
                    term = ''                   
                    while True:
                        term += words[j] + ' '
                        if words[j][-1:] == '"':
                            term = term[:-2]
                            break
                        j += 1
                        assert j < len(words), 'Bad file!'
                    term = '"' + term + '"'
                    terms.append(term)
                    
                #Handle single quotes
                elif words[j][0] == "'":
                    words[j] = words[j][1:]
                    term = ''                   
                    while True:
                        term += words[j] + ' '
                        if words[j][-1:] == "'":
                            term = term[:-2]
                            break
                        j += 1
                        assert j < len(words), 'Bad file!' + str(words)
                    term = "'" + term + "'"
                    terms.append(term)
                else:
                    terms.append(words[j])
                j += 1
            i += 1
        
    #Convert the list of terms into a dictionary
    data = 'data_none'
    alldict = DataModelDict()

    i = 0
    while i < len(terms):
        
        #Single key-value pairs
        if terms[i][0] == '_':
            if len(terms[i+1]) > 1 and terms[i+1][0] == '"' and terms[i+1][-1] == '"':
                terms[i+1] = terms[i+1][1:-1]
            elif len(terms[i+1]) > 1 and terms[i+1][0] == "'" and terms[i+1][-1] == "'":
                terms[i+1] = terms[i+1][1:-1]
            try:
                alldict[data][terms[i]] = terms[i+1]
            except:
                alldict[data] = DataModelDict()
                alldict[data][terms[i]] = terms[i+1]
            i += 2
        
        #Loops
        elif terms[i] == 'loop_':
            i += 1
            if i == len(terms):
                raise ValueError('Bad file!')
            key_list = []
            
            #Retrieve key names
            while terms[i][0] == '_':
                key_list.append(terms[i])
                try:
                    alldict[data][terms[i]] = []
                except:
                    alldict[data] = DataModelDict()
                    alldict[data][terms[i]] = []
                i += 1
                assert i < len(terms), 'Bad file!'
            
            #Save values to key names
            j = 0        
            while i < len(terms) and terms[i] != 'loop_':
                if len(terms[i]) > 0 and terms[i][0] == '_':
                    break 
                if len(terms[i]) > 5 and terms[i][:5] == 'data_':
                    break
                if len(terms[i]) > 1 and terms[i][0] == '"' and terms[i][-1] == '"':
                    terms[i] = terms[i][1:-1]
                elif len(terms[i]) > 1 and terms[i][0] == "'" and terms[i][-1] == "'":
                    terms[i] = terms[i][1:-1]    
                alldict[data][key_list[j]].append(terms[i])
                i += 1
                j += 1
                if j == len(key_list):
                    j = 0
            assert j == 0, 'Bad file! ' + str(key_list)
        
        #Change data block
        elif terms[i][:5] == 'data_':
            data = terms[i]
            i += 1
        
        else:
            print 'Unknown term', terms[i]
            i += 1     
            
    return alldict        

def __calc(symm, site):
    pos = np.empty(3)
    i = 0
    for i in xrange(3):
        terms = []
        c = 0
        while c < len(symm[i]):
            if symm[i][c] == 'x':
                terms.append(site[0])
                c += 1
            elif symm[i][c] == 'y':
                terms.append(site[1])
                c += 1
            elif symm[i][c] == 'z':
                terms.append(site[2])
                c += 1
            elif symm[i][c] == ' ':
                c += 1
            else:
                terms.append(symm[i][c])
                if symm[i][c].isdigit() or symm[i][c] == '.':
                    c += 1
                    while c < len(symm[i]) and (symm[i][c].isdigit() or symm[i][c] == '.'):
                        terms[-1] += symm[i][c]
                        c += 1
                    terms[-1] = float(terms[-1])
                else: 
                    c += 1
        
        #Remove leading signs
        if terms[0] == '+':
            terms = terms[1:]
        elif terms[0] == '-':
            terms[1] = -terms[1]
            terms = terms[1:]
        
        #Multipy and divide
        while '*' in terms or '/' in terms:
            for j in xrange(len(terms)):
                if terms[j] == '*':
                    value = terms[j-1] * terms[j+1]
                    newterms = terms[:j-1] + [value] + terms[j+2:]
                    break
                elif terms[j] == '/':
                    value = terms[j-1] / terms[j+1]                        
                    newterms = terms[:j-1] + [value] + terms[j+2:]
                    break
            terms = newterms 
        
        #add and subtract
        while '+' in terms or '-' in terms:
            for j in xrange(len(terms)):
                if terms[j] == '+':
                    value = terms[j-1] + terms[j+1]
                    newterms = terms[:j-1] + [value] + terms[j+2:]
                    break
                elif terms[j] == '-':
                    value = terms[j-1] - terms[j+1]
                    newterms = terms[:j-1] + [value] + terms[j+2:]
                    break
            terms = newterms

        #Check that all calculations are done
        assert len(terms) == 1, terms            
        
        #Save value and put inside box
        pos[i] = terms[0]          
        while pos[i] >= 1:
            pos[i] -= 1.0
        while pos[i] < 0:
            pos[i] += 1.0  
    return pos

def __tofloat(value):
    """Converts string (possibly with errors) to float."""
    try:
        p_index = value.index('(')
    except:
        p_index = len(value)
    
    return float(value[:p_index])     
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            ucell, elements = cif(f)
        print 'elements =', elements
        print ucell