import math
funlist = ['abs', 'acos', 'atan', 'asin', 'sin', 'tan', 'cos', 'log', 'ln']

def EquationSolver(eqstr, Xvalue, Yvalue=None):
    """(str, number [, number]) -> number

    solves a raw equation string with given Xvalue. Returns number solution.

    >>> EquationSolver('54log(x)', 1)
    0

    >>> EquationSolver('34x+5^(x)(90*8/x)', 3.3)
    44311.92536

    >>> EquationSolver('3xpi', 88.8)
    836.9202829

    >>> EquationSolver('4.4x(2/3)^(2/3)', 90203)
    302886.1992
    """
    #print(Yvalue)
    eqLst = ConvertToList(eqstr, Xvalue, Yvalue)
    solLst = EquationSolver_lo(eqLst)
    return solLst[0]

def ConvertToList(eqstr, Xvalue, Yvalue=None):
    """(str, number [, number]) -> list of items

    Converts an equation string into a list of its individual components,
    replacing 'x' with given Xvalue.
    Parenthesis are listed as string
    Numbers are listed as flost or int
    Operations are listed as string

    >>> ConvertToList('54log(x)', 1)
    [54.0, 'log', '(', 1, ')']

    >>> ConvertToList('3454.56x+5^x(90*8/x)', 3.3)
    [3454.56, 3.3, '+', 5.0, '^', 3.3, '(', 90.0, '*', 8.0, '/', 3.3, ')']

    >>> ConvertToList('3xpi', 88.8)
    [3, 88.8, 3.141592653589793]

    >>> ConvertToList('4.4xe^log(pi)', 90203)
    [4.4, 90203, 2.718281828459045, '^', 'log', '(', 3.14159265
    3589793, ')']

    Examples to test negativehandler:
    >>> ConvertToList('-5sinx', 4.4)
    [-5, 'sin', 4.4]

    >>> ConvertToList('-(5^x)/(-(-5pix))',3)
    [-1, '*', '(', 5, '^', 3, ')', '/', '(', -1, '*', '(', -5, 3.14159265
    3589793, 3, ')', ')']

    >>> ConvertToList('(-5-5x)-(9-x)', 9.0)
    ['(', -5, '-', 5, 9.0, ')', '-', '(', 9, '-', 9.0, ')']

    """

    eqlst = []
    letters = 'qwertuioplkjhgfdsazcvbnm'
    numbers = '1234567890.'
    operators = '!^*-+/'
    parenthesis = '()'
    variables = 'xy'

    numstring = ''
    letterstring = ''

    i = 0
    length = len(eqstr)
    
    while i < length:
##        print(i)
##        print(eqlst)
        if eqstr[i] in letters:
            while i < length and eqstr[i] in letters:
                letterstring += eqstr[i]
                i += 1
            eqlst.append(letterstring)
            letterstring = ''
        
        elif eqstr[i] in numbers:
            while i < length and eqstr[i] in numbers:
                numstring += eqstr[i]
                i += 1
            eqlst.append(float(numstring))
            numstring = ''     

        elif eqstr[i] in operators:
            eqlst.append(eqstr[i])
            i += 1
            
        elif eqstr[i] in variables:
            if eqstr[i] is 'x':
                eqlst.append(Xvalue)
            elif eqstr[i] is 'y':
                eqlst.append(Yvalue)
            i += 1

        elif eqstr[i] in parenthesis:
            eqlst.append(eqstr[i])
            i += 1
            
    for i in range(0, len(eqlst)):
        if eqlst[i] == 'pi':
            eqlst[i] = math.pi
        elif eqlst[i] == 'e':
            eqlst[i] = math.e
    
    eqlst = negativehandler(eqlst)
    #print(eqlst)
    return eqlst
                       
def negativehandler(eqlst):
    """(list of str) -> list of str
    This function handles negative signs
    
    If the operator is a '-', it must be diferenciated as
    an operator or a negative sign.
    
    WARNING: A negative sign must either be at the front of the equation
    or be preceded by an open parenthesis.
    """
    
    returnlst = []

    letters = 'qwertyuioplkjhgfdsazcvbnm'
    parenthesis = '()'
    
    i=0
    length = len(eqlst)
    
    while i < length:
        item = eqlst[i]
    
        if str(item) not in '-':
            returnlst.append(item)
            i += 1
            
        elif str(item) in '-':
            
            if i == 0 or eqlst[i-1] == '(':
                obj = eqlst[i+1]
                
                #if the negative sign is followed by a number
                if type(obj) is int or type(obj) is float:
                    
                    Negative_num = -1*obj
                    returnlst.append(Negative_num)
                    i += 2
                #if the negative sign is followed by a parenthesis
                elif obj in parenthesis:
                    returnlst.append(-1)
                    returnlst.append('*')
                    returnlst.append('(')
                    i += 2                 
                #if the negative sign is followed by a special letter function.
                elif obj in funlist:
                    
                    returnlst.append(-1)
                    returnlst.append('*')
                    returnlst.append(obj)
                    i += 2
                             
            else:
                #not a negative sign
                returnlst.append(item)
                i += 1
        
    return returnlst
        
def EquationSolver_lo(eqLst):
    """(list of str) -> float

    Solves an equation with individual parts split into a list of str.
    all parenthesis must be string
    all operations must be string.
    all numbers must be int or float.

    >>> list = ['(', '(', 3, ')', ')']
    >>> EquationSolver(list)
    [3]

    >>> list = [ '(', 3, '+', math.pi, ')', '/', '(',
    2, '^', '(', 4, '+', 1, ')', ')']
    >>> EquationSolver(list)
    [0.19192477042468103]

    >>> list = ['ln', '(', math.e, '^', '(', 5, '+', 1, '^',
    '(', 'sin', '(', math.pi, ')', ')', ')', ')']
    >>> EquationSolver(list)
    [6.0]

    """
    
    parstring = ''
    
    while parSearch(eqLst):
        #print('while')
        for i in range(0, len(eqLst)):
            if str(eqLst[i]) not in '()':
                i += 1
                #print('if')
            elif str(eqLst[i]) in '()':
                #print('elif')
                if parstring == '':
                    parstring = eqLst[i]
                    parindex = i
                    i += 1
                elif parstring == eqLst[i]:
                    parindex = i
                    i += 1
                elif parstring != eqLst[i]:
                    solution = EquationSolver_np(eqLst[parindex+1:i])
                    #print(parindex)
                    #print(solution)
                    eqLst = lstReplace(eqLst, parindex, i+1, solution)
                    #print(eqLst)
                    parstring = ''
                    break
    if len(eqLst) == 1:
        return eqLst

    solution = EquationSolver_np(eqLst)
    return solution
                    

def EquationSolver_np(eqLst):
    """(list of str) -> float

    Solves an equation with individual parts split into a list of str.
    all operations must be string.
    all numbers must be int or float.
    Parenthesis must be resolved.

    >>> list = [3]
    >>> EquationSolver_np(list)
    [3]

    >>> list = [3,'/',2]
    >>> EquationSolver_np(list)
    [1.5]

    >>> list = [4, '+', 5.5]
    >>> EquationSolver_np(list)
    [9.5]
    
    >>> list = ['sin', math.pi, 5, 6,]
    >>> EquationSolver_np(list)
    [0]

    >>> list = [5, '^', 6, 5, '/', 5, '+', 'tan', math.pi, math.pi, 6]
    >>> EquationSolver_np(list)
    [15625]

    >>> list = [5, '+', 5, 5, 5, '-', 5, '/',
    5, '/', 5, 5, '/', 5, '+', 5, '-', 5, '+', 5, 5, 5, '/', 5, '/', 5]
    >>> EquationSolver_np(list)
    [134.8]

    """    

    if len(eqLst) == 1:
        return eqLst

    #Resolve special math functions.
    solutionlist1 = []
    i = 0
    length = len(eqLst)
    while i < length:
        if eqLst[i] in funlist:
            fun = eqLst[i]
            var = eqLst[i+1]
            if fun == 'abs':
                value = abs(var)
            elif fun == 'acos':
                value = math.acos(var)
            elif fun == 'asin':
                value = math.asin(var)
            elif fun == 'atan':
                value = math.atan(var)
            elif fun == 'sin':
                value = math.sin(var)
            elif fun == 'cos':
                value = math.cos(var)
            elif fun == 'tan':
                value = math.tan(var)
            elif fun == 'log':
                value = math.log10(var)
            elif fun == 'ln':
                value = math.log(var)
            solutionlist1.append(value)
            i += 2
        else:
            solutionlist1.append(eqLst[i])
            i += 1

    #print(solutionlist1)
    if len(solutionlist1) == 1:
        return solutionlist1

    #Resolve exponents        
    solutionlist2 = []
    length = len(solutionlist1)
    i = 0
    while i < length:
        if i == length - 1:
            solutionlist2.append(solutionlist1[i])
            break
        if eqLst[i] == '^':
            value = math.pow(eqLst[i-1], eqLst[i+1])
            solutionlist2.append(value)
            i += 2
        elif eqLst[i+1] == '^':
            i += 1
        else:
            solutionlist2.append(solutionlist1[i])
            i += 1
    
    #print(solutionlist2)
    if len(solutionlist2) == 1:
        return solutionlist2

    #Insert '*' where nultiplication is understood.
    solutionlist3 = []
    length = len(solutionlist2)
    i = 0
    while i < length:
        if i < length-1 and str(solutionlist2[i]) not in '!^*-+/' and str(solutionlist2[i+1]) not in '!^*-+/':
            j = i
            while i < length and str(solutionlist2[i]) not in '!^*-+/':
                i += 1
            for index in range(j,i-1):
                solutionlist3.append(solutionlist2[index])
                solutionlist3.append('*')
            solutionlist3.append(solutionlist2[i-1])
        
        else:
            solutionlist3.append(solutionlist2[i])
            i += 1

    #print(solutionlist3)
    if len(solutionlist3) == 1:
        return solutionlist3

    #Resolve multiplication and division.    
    solutionlist4 = []
    length = len(solutionlist3)
    i = 0
    while i < length:
        
        if i < length-1 and str(solutionlist3[i]) in '*/':
            j = i
            fun = solutionlist3[i]
            if fun in '/':
                value = solutionlist3[i-1]/solutionlist3[i+1]
            elif fun in '*':
                value = solutionlist3[i-1]*solutionlist3[i+1]
            
            while i < length and solutionlist3[i] in '*/':
                i += 2
            if j + 2 == i:
                solutionlist4.append(value)
                continue
            for index in range(j + 2,i,2):
                fun = solutionlist3[index]
                if fun in '/':
                    value = value/solutionlist3[index+1]
                elif fun in '*':
                    value = value*solutionlist3[index+1]
                    
            solutionlist4.append(value)
            
##        elif solutionlist3[i] in '/':
##            value = solutionlist3[i-1]/solutionlist3[i+1]
##            solutionlist4.append(value)
##            i += 2
##        elif solutionlist3[i+1] in '*/':
##            i += 1
        elif i < length-1 and str(solutionlist3[i+1]) in '*/':
            i += 1 
            
        else:
            solutionlist4.append(solutionlist3[i])
            i += 1

    #Handle negatives
##    if solutionlist4[0] == '-':
##        solutionlist4 = negativehandler(solutionlist4)

    #print(solutionlist4)
    if len(solutionlist4) == 1:
        return solutionlist4

    #Resolve addition and subtraction.
    solutionlist = []

    fun = solutionlist4[1]

    if fun in '+':
        #print(solutionlist4[0], '+', solutionlist4[2])
        value = solutionlist4[0] + solutionlist4[2]
    elif fun in '-':
        value = solutionlist4[0] - solutionlist4[2]
    if len(solutionlist4) == 3:
        solutionlist.append(value)
        return solutionlist
    else:
        for index in range(3,len(solutionlist4), 2):
            fun = solutionlist4[index]
            if fun in '+':
                value = value + solutionlist4[index+1]
            elif fun in '-':
                value = value - solutionlist4[index+1]
    solutionlist.append(value)

    return solutionlist
             
def lstReplace(lst, iStart, iEnd, replacemt):
    """(list of items, int, item) -> list of items
    Returns a corrected version of lst with
    lst[iStart:iEnd] replaced with replacemt.

    >>> list = [0,1,2,3,4,5]
    >>> lstReplace(list,3,5, 'x')
    [0,1,2,'x',5]

    """
    
    returnList = lst[:iStart]
    returnList += replacemt
    if iEnd == len(lst):
        return returnList
    
    returnList += lst[iEnd:]

    return returnList
    
def parSearch(lst):
    """(list of items) -> bool

    Indicates whether there is a parenthesis in lst.

    """

    for item in lst:
        if str(item) in '()':
            return True
        
    return False
