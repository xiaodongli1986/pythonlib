####n##########################
### Arrays
###############################

# print
def print_2dmat(X, fmt='%.5f', div = ' '):
    for x in X:
        nowstr = ''
        for xx in x:
            nowstr += fmt%xx
            nowstr += div
        print( nowstr)

# Transpose 2d array
def transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]
# Reverse 1d array
def reverse_1darray(A):
    return [A[i] for i in range(len(A)-1,-1,-1)]

# Some manipulation
def XmultY(X, Y, a=1,b=1,c=1):
        ''' (X[row]**a)*(Y[row]**b)*c for row in range(len(X)) '''
        if a==1 and b==1 and c==1:
                return [X[row]*Y[row] for row in range(len(X))]
        else:
                return [(X[row]**a)*(Y[row]**b)*c for row in range(len(X))]

def Xmulty(X, y, a=1,b=1,c=1):
        ''' (X[row]**a)*(y**b)*c for row in range(len(X)) '''
        if a==1 and b==1 and c==1:
                return [X[row]*y for row in range(len(X))]
        else:
                return [(X[row]**a)*(y**b)*c for row in range(len(X))]


def Xdivy(X, y):
        ''' (X[row]**a)/y for row in range(len(X)) '''
        return [X[row]/y for row in range(len(X))]



def XplusY(X, Y, a=1,b=1):
        if a==1 and b==1:
                return [X[row]+Y[row] for row in range(len(X))]
        else:
                return [(X[row]*a)+(Y[row]*b) for row in range(len(X))]

def product(iterable):
    return reduce(operator.mul, iterable, 1)

def sumarray(X,dim=1):
        rlt = 0
        if dim == 1:
                rlt = sum(X)
        elif dim == 2:
                rlt = sum([sum(xx) for xx in X])
        elif dim == 3:
                rlt = sum([sum([sum(xx) for xx in xxx]) for xxx in X])
        elif dim == 4:
                rlt = sum([sum([sum([sum(xx) for xx in xxx]) for xxx in xxxx]) for xxxx in X])
        elif dim == 5:
                rlt = sum([sum([sum([sum([sum(xx) for xx in xxx]) for xxx in xxxx]) for xxxx in xxxxx]) for xxxxx in X])
        return rlt

def subarray_4from5(X,idim=0):
        len1, len2, len3, len4, len5 = len(X), len(X[0]), len(X[0][0]), len(X[0][0][0]), len(X[0][0][0][0])
        if idim ==0:
                return [[[[[X[row][row1][row2][row3][row4] for row1 in range(len2)] for row2 in range(len3)] for row3 in range(len4)] for row4 in range(len5)] for row in range(len1)]
        elif idim ==1:
                return [[[[[X[row1][row][row2][row3][row4] for row1 in range(len1)] for row2 in range(len3)] for row3 in range(len4)] for row4 in range(len5)] for row in range(len2)]
        elif idim ==2:
                return [[[[[X[row1][row2][row][row3][row4] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len4)] for row4 in range(len5)] for row in range(len3)]
        elif idim ==3:
                return [[[[[X[row1][row2][row3][row][row4] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len3)] for row4 in range(len5)] for row in range(len4)]
        elif idim ==4:
                return [[[[[X[row1][row2][row3][row4][row] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len3)] for row4 in range(len4)] for row in range(len5)]

def sum4of5darray(X,idim=0):
        len1, len2, len3, len4, len5 = len(X), len(X[0]), len(X[0][0]), len(X[0][0][0]), len(X[0][0][0][0])
        #print
        #print len1, len2, len3, len4, len5
        if idim ==0:
                return [sumarray(XX,4) for XX in X]
        elif idim ==1:
                return [sumarray([[[[X[row1][row][row2][row3][row4] for row1 in range(len1)] for row2 in range(len3)] for row3 in range(len4)] for row4 in range(len5)],4) for row in range(len2)]
        elif idim ==2:
                return [sumarray([[[[X[row1][row2][row][row3][row4] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len4)] for row4 in range(len5)],4) for row in range(len3)]
                #return [sumarray(X[0:len1][0:len2][row][0:len4][0:len5],4) for row in range(len3)]
        elif idim ==3:
                return [sumarray([[[[X[row1][row2][row3][row][row4] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len3)] for row4 in range(len5)],4) for row in range(len4)]
                #return [sumarray(X[0:len1][0:len2][0:len3][row][0:len5],4) for row in range(len4)]
        elif idim ==4:
                return [sumarray([[[[X[row1][row2][row3][row4][row] for row1 in range(len1)] for row2 in range(len2)] for row3 in range(len3)] for row4 in range(len4)],4) for row in range(len5)]
                #return [sumarray(X[0:len1][0:len2][0:len3][0:len4][row],4) for row in range(len5)]
# Selecting columns from multi dim arrau
def Xfromdata(data, i=0):
    return [data[row][i] for row in range(len(data))]
def XYfromdata(data, i=0, j=1):
    return [ data[row][i] for row in range(len(data))], [ data[row][j] for row in range(len(data))]    
def XYZfromdata(data, i=0, j=1, k=2):
    return [ data[row][i] for row in range(len(data))], [ data[row][j] for row in range(len(data))],\
        [ data[row][k] for row in range(len(data))] 
def Xsfromdata(data, ilist):
    return ([ data[row][ilist[i]] for row in range(len(data))] for i in range(len(ilist)))

def Xsfrom2ddata(data, ilist):
    return ([[ data[row1][row2][ilist[i]] for row2 in range(len(data[row1]))] for row1 in range(len(data))] for i in range(len(ilist)))

def Xfrom2ddata(data, i=0):
    return [[ data[row1][row2][i] for row2 in range(len(data[row1]))] for row1 in range(len(data))]


    
####
def meanY(Y): ### mean value of an array
    return sum(Y)/float(len(Y))
def meanPartY(Y,i1,i2): ### mean value of an array
    return sum(Y[i1:i2+1])/float((i2-i1+1))

def randpart(A, randrat = 1.0):
    B = []
    for X in A:
        if random.uniform(0,1) < randrat:
            B.append([x for x in X])
    return B

# Something about index; random selection from an array;...
def Xs(X, index, deepcopy=False):
        if not deepcopy:
                return [X[row] for row in index]
        else:
                return [copy.deepcopy(X[row]) for row in index]

def randrange(imax, imin=0, rat=0.1):
        rlt = []
        for i in range(imin, imax):
                if random.uniform(0,1) < rat:
                        rlt.append(i)
        return rlt

def randXpart(X, rat=0.1, deepcopy=False, onlyreturnvalue=False, onlyreturnindex=False):
        rltindex = randrange(len(X), 0, rat = rat)
        rltX = Xs(X, rltindex, deepcopy=deepcopy)
        if onlyreturnvalue:
                return rltX
        elif onlyreturnindex:
                return rltindex
        else:
                return rltX, rltindex

### averaging of 2d arrays... arraylists will be 3d list; first index is the index of array.
def get_avg_array2D(arraylists, furtherdivfac = 1.0):
    import copy
    avgarray = copy.deepcopy(arraylists[0])
    divfac   = float(len(arraylists)) * furtherdivfac
    nrow = len(arraylists[0])
    for i in range(nrow):
        ncol = len(arraylists[0][i])
        for j in range(ncol):
            for iarray in range(1,len(arraylists)):
                avgarray[i][j] += arraylists[iarray][i][j]
            avgarray[i][j] /= divfac
    return avgarray
def get_std_array2D(arraylists, furtherdivfac = 1.0):
    import copy
    return [[ get_stat_from_list([ arraylists[row][row1][row2]  for row in range(len(arraylists))])[2] \
                for row2 in range(len(arraylists[0][0]))  ] for row1 in range(len(arraylists[0]))]        
def get_avg_array(arraylists, furtherdivfac = 1.0):
    import copy
    avgarray = copy.deepcopy(arraylists[0])
    divfac   = float(len(arraylists)) * furtherdivfac
    nrow = len(arraylists[0])
    for i in range(nrow):
        for iarray in range(1,len(arraylists)):
            avgarray[i] += arraylists[iarray][i]
        avgarray[i] /= divfac
    return avgarray
def get_avgstd_array(arraylists, ConsiderVarVar=False):
    nrow = len(arraylists[0])
    avgarray, avgstdarray = [0 for row in range(nrow)], [0 for row in range(nrow)]
    for i in range(nrow):
        X = []
        for iarray in range(len(arraylists)):
            X.append(arraylists[iarray][i])
        if not ConsiderVarVar:
                av, var, aver = get_stat_from_list(X)
        else:
                av, var, aver, varer = get_stat_from_list(X,getvarer = True)
                var = var+varer
                aver = np.sqrt(var) / np.sqrt(len(X)-1.0)
        avgarray[i] = av; avgstdarray[i] = aver*np.sqrt(float(len(arraylists))-1.0)
    return avgarray, avgstdarray

### conversion between 1d/2d arrays
def get_1darray_from_2d(array_2d):
        array_1d=[]
        for row in range(len(array_2d)):
                array_1d += array_2d[row]
        return array_1d
    
def get_2darray_from_1d(A, numi, numj):
    B = [[0 for row2 in range(numj)] for row1 in range(numi)]
    k = 0
    for i in range(numi):
        for j in range(numj):
            B[i][j] = A[k]
            k+=1
    return B

# Somthing about select range
def Xinrange(X, xmin, xmax, onlyreturnvalue=False, onlyreturnindex=False):
        rltX, rltindex = [], []
        for row in range(len(X)):
                if xmin <X[row] < xmax:
                        rltX.append(X[row]); rltindex.append(row)
        if onlyreturnvalue:
                return rltX
        elif onlyreturnindex:
                return rltindex
        else:
                return rltX, rltindex
def list_Xinrange(X, xlist, onlyreturnvalue=False, onlyreturnindex=False):
        rlt = []
        for row in range(len(xlist)-1):
                rlt.append(Xinrange(X, xlist[row], xlist[row+1], onlyreturnvalue=onlyreturnvalue, onlyreturnindex=onlyreturnindex))
        return rlt
def XYinrange(X, Y, xmin, xmax, ymin, ymax, onlyreturnvalue=False, onlyreturnindex=False):
        rltX, rltY, rltindex = [], [], []
        for row in range(len(X)):
                if xmin <X[row] < xmax and ymin < Y[row] < ymax:
                        rltX.append(X[row]); rltindex.append(row); rltY.append(Y[row]); 
        if onlyreturnvalue:
                return rltX, rltY
        elif onlyreturnindex:
                return rltindex
        else:
                return rltX, rltY, rltindex
        

# Something about Func Selection...
def FuncSelcX(X, func):
        rlt = []
        for x in X:
                if func(x):
                        rlt.append(x)
        return rlt
def FuncSelcX_index(X, func):
        rlt = []
        for row in range(len(X)):
                if func(X[row]):
                        rlt.append(row)
        return rlt

# Polynomial fitting
        
def polyfitY(X, Y, deg):
    polyfitrlt = polyfit(X, Y, deg=deg)
    return polyval(polyfitrlt, X) 
      
def polystr(polyfit, deg=2, fmt='%.1f', valstr = 'x'):
            nowstr = ''
            nowdeg = deg
            for i in range(deg):
                coef = polyfit[i]
                if i>= 1 and coef > 0:
                    nowstr += ' + '
                else:
                    nowstr += '  '
                nowstr += (fmt%coef)
                if nowdeg >1:
                    nowstr += (' '+valstr+'^'+str(nowdeg))
                else:
                    nowstr += (' '+valstr)
                nowdeg -= 1
            i = deg
            coef = polyfit[i]
            if i>= 1 and coef > 0:
                    nowstr += ' + '
            nowstr += (fmt%coef)
            return nowstr


### return rows where the redshift lies within certain region
### universially can be applied to any multi-component array to limit the range of a certain component
def zlimited_data(bossdata, zmin, zmax, zcol=2):
    newdata = Xs(bossdata, Xinrange(Xfromdata(bossdata, zcol), zmin, zmax, onlyreturnindex=True))
    return newdata

### sort a multi-component array; adopting to a certain component as the key
def sortXs(Xs, col=0):
    def sortedkey(X,iptcol=col):
        return X[iptcol]
    return sorted(Xs, key=sortedkey)

def id_in_array(x, X, tol=0.001): ## id in the array
    nowrow = -1
    for row in range(len(X)):
        if abs(x-X[row])<tol:
            nowrow = row; return nowrow
    return nowrow
def ArraytoStr(Array, div = ' '): ## convert array to string
    str0 = ''
    for x in Array:
        str0 += (str(x)+div)
    return str0
def sumlist(A):
        rlt = []
        for a in A:
                rlt += a
        return rlt
def array_merge(A):
        return sumlist(A)
def get_absarray(X):
        return [abs(x) for x in X]
def get_sumarray(X,Y):
        return [Y[row]+X[row] for row in range(len(Y))]
def get_diffarray(X,Y):
        return [Y[row]-X[row] for row in range(len(Y))]
def get_divarray(X,Y):
        return [X[row]/Y[row] for row in range(len(Y))]

def get_divarray_2d(X,Y):
        return [[X[row1][row2]/Y[row1][row2] for row2 in range(len(Y[row1]))] for row1 in range(len(Y))]

def get_diffarray_3darray(X,Y,misplace=None):
        if misplace==None:
                return [[[ Y[row1][row2][row3] - X[row1][row2][row3] for row3 in range(len(X[row1][row2]))] \
                        for row2 in range(len(X[row1]))] for row1 in range(len(X))]
        else:
                return [[[ Y[row1][row2][row3] - X[row1-misplace][row2][row3] for row3 in range(len(X[row1][row2]))] \
                        for row2 in range(len(X[row1]))] for row1 in range(len(X))]
### functions that normalize an array
# averaged value of X
def avgarray(X,wei=[]):
        if wei == []:
                return sum(X)/float(len(X))
        else:
                return sum([X[row]*wei[row] for row in range(len(X))]) / float(sum(wei))
def avgarray_2d(X,wei=[]):
        if wei == []:
                return [sum([X[i][row] for i in range(len(X))]) / float(len(X)) for row in range(len(X[0]))]
        else:
                sumwei = sum(wei)
                return [sum([X[i][row]*wei[i] for i in range(len(X))]) / sumwei for row in range(len(X[0]))]

def shiftozero(X,wei=[],returnavg=False):
        avg = avgarray(X,wei);
        if not returnavg:
                return [X[row]-avg for row in range(len(X))]
        else:
                return avg, [X[row]-avg for row in range(len(X))]
def normto1(X,wei=[],returnavg=False):
        avg = avgarray(X,wei);
        if not returnavg:
                return [X[row]/avg for row in range(len(X))]
        else:
                return avg, [X[row]/avg for row in range(len(X))]
def normtolike(X, xmin=None):
        ''' Normalize a chisq-array to a likelihood array
        if xmin == None: xmin = min(X)
        return [exp(-0.5*(xx-xmin) for xx in X]'''
        if xmin == None: xmin = min(X)
        return [exp(-0.5*(xx-xmin)) for xx in X]
def norm2dcovmat(A):
    ''' normalize covariance matrix so that diagnal elements equal to 1'''
    nA = len(A)
    diag = [A[row][row] for row in range(nA)]
    B = [[A[i][j]/np.sqrt(diag[i]*diag[j]) for j in range(nA)] for i in range(nA)]
    return B
def normto1_2darray(X,Skip_element=[],wei=[],returnavg=False):
        '''
        Skip_element = [str(xx[0])+'_'+str(xx[1]) for xx in Skip_element]
        sumX = 0; numX = 0
        for row1 in range(len(X)):
                for row2  in range(len(X[0])):
                        nowrowstr = str(row1)+'_'+str(row2)
                        if nowrowstr in Skip_element: continue
                        sumX  += X[row1][row2]
                        numX +=1
        numX += 0.0
        avgX = sumX / numX
        normedX = [[xxx/avgX for xxx in xx] for xx in X ]
        if not returnavg:
                return normedX
        else:
                return avgX, normedX'''
        Skip_element = [str(xx[0])+'_'+str(xx[1]) for xx in Skip_element]
        sumX = 0; numX = 0
        for row1 in range(len(X)):
                for row2  in range(len(X[0])):
                        nowrowstr = str(row1)+'_'+str(row2)
                        if nowrowstr in Skip_element: continue
                        sumX  += X[row1][row2]
                        numX +=1
        numX += 0.0
        avgX = sumX / numX
        normedX = [[xxx/avgX for xxx in xx] for xx in X ]
        if not returnavg:
                return normedX
        else:
                return avgX, normedX
def norm_power(X, deg=0.9):
                        avg = sum([x**deg for x in X])**(1.0/deg) / float(len(X))
                        return [x/avg for x in X]
def norm_OneSkip(X, deg=1.0):
                        avg = sum(X) / float(len(X))
                        imid = len(X)/2 + 1
                        Y = [X[row]/avg for row in range(imid)] + [X[row]/avg for row in range(imid+1,len(X))]
                        return Y
def norm_SkipLast(X):
    return normto1(X)[0:len(X)-1]
def norm_SkipFirst(X):
    return normto1(X)[1:len(X)]
def norm_powerC(X, deg=1.0):
                        avg = sum([x**deg for x in X])**(1.0/deg) / float(len(X))
                        return [X[row]/avg for row in range(len(X)-2)]



def arraysq(X):
        return [x**2.0 for x in X]

####
## output/loadin 3d array (each row is 2d)
def output_3darray_to_file_eachrowis2d(covmats, filename):
    nowf = open(filename, 'w')
    numrbin = len(covmats)
    for irbin in range(numrbin):
        nowstr = array_to_str(get_1darray_from_2d(covmats[irbin])) + '\n'
        nowf.write(nowstr)

def get_3darray_from_file_eachrowis2d(filename, numi, numj):
    data = np.loadtxt(filename)
    return [get_2darray_from_1d(data[row], numi, numj) for row in range(len(data))]



### Search for the point, above/below which the sum of elements firstly exceed the goal; 

def find_accumulated_point(A, numgoal, findtype = 'G'):#G means greater; L means less
    na = len(A)
    if findtype == 'G':
        ip = na-1; num = A[na-1]
        while num < numgoal:
            ip -= 1
            num += A[ip]
        #realip = ip + 1.0/(float(A[ip])) * (numgoal - num + A[ip])        
    if findtype == 'L':
        ip = 0; num = A[0]
        while num < numgoal:
            ip += 1
            num += A[ip]
        #realip = ip + 1.0/(float(A[ip])) * (numgoal - num + A[ip])
    #print "num - A[ip], num, numgoal, realip", num - A[ip], num, numgoal, realip
    return ip, num

###############################
## Conversion between numbers and characters/strings;
###############################
### convert an array to str

def str_to_numbers(str1,exitcode='#', do_float_conver=True):
    floatlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        if duru == exitcode:
            break
        elif duru == '\n' or duru == ' ' or duru == '\t' or duru == ',':
            if numberunderconstruct == True:
                floatlist.append(str2); str2 = ''
                numberunderconstruct = False
        else:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        floatlist.append(str2);
    if do_float_conver:
        floatlist = [float(floatlist[row]) for row in range(len(floatlist))]
    return floatlist

def get_str_list(str1, div = [' ', ',', '\t', '\n']):
    strlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        strisdiv = False
        for sepstr in div:
            if duru == sepstr:
                strisdiv = True
                if numberunderconstruct == True:
                    strlist.append(str2); str2 = ''
                    numberunderconstruct = False
        if strisdiv == False:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        strlist.append(str2);
    return strlist

def str_to_numbers(str1,exitcode='#', do_float_conver=True):
    floatlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        if duru == exitcode:
            break
        elif duru == '\n' or duru == ' ' or duru == '\t' or duru == ',':
            if numberunderconstruct == True:
                floatlist.append(str2); str2 = ''
                numberunderconstruct = False
        else:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        floatlist.append(str2);
    if do_float_conver:
        floatlist = [float(floatlist[row]) for row in range(len(floatlist))]
    return floatlist

def get_str_list(str1, div = [' ', ',', '\t', '\n']):
    strlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        strisdiv = False
        for sepstr in div:
            if duru == sepstr:
                strisdiv = True
                if numberunderconstruct == True:
                    strlist.append(str2); str2 = ''
                    numberunderconstruct = False
        if strisdiv == False:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        strlist.append(str2);
    return strlist

def array_to_str(A,div=' '):
    nowstr = '';
    for a in A:
        nowstr += str(a)
        nowstr += div
    return nowstr

def array_to_str_fmtted(A,div=' ', fmt='%.4f'):
    nowstr = '';
    for a in A:
        nowstr += fmt%a
        nowstr += div
    return nowstr

def strarray_to_float1D(A, replace = False):
        if not replace:
                return [float(a) for a in A]
        else:
                for row in range(len(A)):
                        A[row] = float(A[row])
                return A

def strarray_to_float2D(AA):
                return [[float(        a) for a in A] for A in AA]



def index_of_max(X):
            im = 0; xmax = X[0];
            for row in range(len(X)):
                if X[row] > xmax:
                    xmax = X[row]; im = row
            return im

def convert_loaded2Ddata_to_array(data):
        return [[data[row1][row2] for row2 in range(len(data[0]))] for row1 in range(len(data))]
        
def get_mid_array1d(data, deg=1):
        if deg == 1:
                X = [(data[row]+data[row+1])*0.5 for row in range(len(data)-1)]
        else:
                X = [((data[row]**deg+data[row+1]**deg)*0.5)**(1.0/deg) for row in range(len(data)-1)]
        return X
        
def get_mid_array2d(data):
        return [[(data[row1][row2]+data[row1][row2+1]+data[row1+1][row2]+data[row1+1][row2+1])*0.25 for row2 in range(len(data[row1])-1)] for row1 in range(len(data)-1)]

# merge two list, within finite digits
def merge_two_list_finitedigits(A,B, numdigit=7,sortlist=True):
    X = [round(x,numdigit) for x in list(A)] + [round(x,7) for x in list(B)]
    X = list(set(X))
    if sortlist: X.sort()
    return X

def check_interval(A):
    return min(abs(A[i]-A[i-1]) for i in range(len(A)))

def array_fracbin(A, nbin=10, i1=None, i2=None):
        if i1==None: i1=0
        if i2==None: i2=len(A)
        edges = np.linspace(i1,i2,nbin+1)
        #print edges
        B = range(nbin)
        #print B
        for i in range(nbin):
                B[i] = 0
                iA, iD = edges[i], edges[i+1]
                iB, iC = int(iA)+1, int(iD)
                if iC >= iB:
                    B[i] = sum(A[iB:iC])
                    if iB>0: B[i] += A[iB-1]*(iB-iA)
                    if iC<len(A): B[i] += A[iC]*(iD-iC)
                elif iC < iB:
                    B[i] = A[iC]*(iD-iA)
        return B
def array_intbin(A, nbin=10, i1=None, i2=None):
        '''Split an array into nbin bins; keep int bins (will not split one element and assign to different bins)'''
        if i1==None: i1=0
        if i2==None: i2=len(A)
        edges = np.linspace(i1,i2,nbin+1)
        edges = [int(edge+0.5) for edge in edges]
        #print edges
        B = range(nbin)
        for i in range(nbin):
                B[i] = 0
                iA, iB = edges[i], edges[i+1]
                #print A[iA:iB]; 
                B[i] = sum(A[iA:iB])
        #print sum(A), sum(B)
        return B
        ''' For test
        A = range(1000)
        nbin = 100
        for i1 in range(2000):
            B = array_intbin(A, nbin, i1)        
            if abs(sum(A[i1:]) - sum(B)) != 0:
                print i1, sum(A[i1:]), sum(B)
                break'''


def drop_duplicate(duplicate):
    final_list = []
    for num in duplicate:
        if num not in final_list:
            final_list.append(num)
    return final_list
