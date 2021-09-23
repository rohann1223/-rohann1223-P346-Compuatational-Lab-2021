#library of functions for computation lab

#reading and writing a matrix from given coefficients

def matrix(r,c): #defining the matrix, asking for rows (r) and coloumns (c) from the user
    mat = []
    for i in range (r):
        row = []
        for j in range (c):
            row.append(float(input ("Enter coefficients of rows ")))  #asking input for rows
        mat.append(row)
        print ("Now write the next row ")  #ending the row and going to the next
    return mat

    for k in mat:
        print (k)

def transpose(T): 
    result =[[0 for c in range(len(T[0]))] for r in range(len(T))]
    for i in range(len(T)):
        # iterate through columns
        for j in range(len(T[0])):
            result[j][i] = T[i][j]
    return result


#partial pivoting

def PartialPivot(Ab, m, rows, cols):
    global n,swapnumber          #global variable to store how many swap are done
    n = 0
    swapnumber = 0
    pivot = Ab[int(m)][int(m)]          #starting pivot of matrix
    for i in range (int(rows)):         
        if pivot < Ab[int(i)][int(m)]:  #checking with other elements of the same coloumn
            pivot = Ab[int(i)][int(m)]
            n += 1
            swapnumber = i
    if swapnumber != 0:
        swapRows(Ab, m, swapnumber, cols)    #swapping if condition satisfies
            
    if int(pivot) == 0:
        print ("No unique solution")   #if pivot is 0 at end it returns no solution
        return None
    
#code for swapping rows
def swapRows(Ab, old, new_r, cols):
    temp = []       #temp list to store old list

    for c in range (0, int(cols)):
        temp.append(Ab[int(old)][c])
        Ab[int(old)][c] = Ab[int(new_r)][c]     #swapping values
        Ab[int(new_r)][c] = temp[c]  

#beginning gauss jordan elimination
def gaussJordan(Ab, rows, cols):
    det = 1                 #to calculate determinant
    PartialPivot(Ab, 0, rows, cols)
    for i in range (rows):
        if Ab == None:    #if pivot = 0
            return Ab
        else:
            fact = Ab[i][i]
            det = det*float(fact)       #calculating determinant by multiplying diagonal elements
            for c in range(i, cols):
                if Ab[i][i] == 0:
                    raise ZeroDivisionError   #if divided by 0
                else:
                    Ab[i][c] = float(Ab[i][c])*1/float(fact)     #reciprocal of the pivot element

            for r in range(rows):
                if (r == i or Ab[r][i]== 0): continue    #if 0 reached or diagonal element
                factor = float(Ab[r][i])
                for j in range(i,cols):
                    Ab[r][j] = float(Ab[r][j])- (factor*float(Ab[i][j]))  #subtracting to obtain 0

    for result in range(rows):
        print (round((Ab[result][cols-1])))   #for printing the result
    return int(det)*((-1)**(n))          #for printing the determinant

    
#code to call function for determinant
def determinant(Ab, rows, cols):
    d = gaussJordan(Ab, rows, cols)
    print ("The determinant is " + str(d))

def inverse(a, r, c):
    x = len(a) #defining the range through which loops will run
    #constructing the n X 2n augmented matrix
    P = [[0.0 for i in range(len(a))] for j in range(len(a))]
    for i in range(3):
        for j in range(3):
            P[j][j] = 1.0
    for i in range(len(a)):
        a[i].extend(P[i])
    #main loop for gaussian elimination begins here
    for k in range(x):
        if abs(float(a[k][k])) < 1.0e-12:
            for i in range(k+1, x):
                if abs(float(a[i][k])) > abs(float(a[k][k])):
                    for j in range(k, 2*x):
                        a[k][j], a[i][j] = a[i][j], a[k][j] #swapping of rows
                    break
        pivot = a[k][k] #defining the pivot
        if pivot == 0: #checking if matrix is invertible
            print("This matrix is not invertible.")
            return
        else:
            for j in range(k, 2*x): #index of columns of the pivot row
                a[k][j] = float(a[k][j])/float(pivot)
            for i in range(x): #index the subtracted rows
                if i == k or a[i][k] == 0: continue
                factor = float(a[i][k])
                for j in range(k, 2*x): #index the columns for subtraction
                    a[i][j] -= factor * a[k][j]
    for i in range(len(a)): #displaying the matrix
        for j in range(n, len(a[0])):
            print("{:.2f}".format(a[i][j]), end = " ") #printing upto 2 places in decimal. 
        print()

#multiplication of two 3x3 matrices
def multiplysquare(x,y):
    product = [[0 for c in range(len((x)[0]))] for r in range (len(x))]
    for i in range(len(x)):
    
        # iterating by column by B
        for j in range(len(x[i])):
    
            # iterating by rows of B
            for k in range(len(x[i])):
                product[i][j] += x[i][k] * y[k][j]
            product[i][j] = round(product[i][j], 6)
    for i in range(len(x)):
        for j in range(len(x[i])):
            if product[i][j] < 10**(-15):
                product[i][j] = 0

    print("The product of ")        #for showing the products
    for a in x:
        print(a)

    print ("and ")

    for b in y:
        print(b)

    print ("is")

    for ab in product:
        print(ab)       #final product
    return product

#writing a matrix using a txt file
def matrixtxt(dir):
    with open(dir, 'r') as f:
        l = [[float(num) for num in line.split(',')] for line in f]
    return l
        
#doolittle algorith to calculate LU
def LUDoolite(M):
    PartialPivot(M, 0, len(M), len(M[0]))  #partial pivoting given matrix
    n = len(M)        
    lower = [[0 for x in range(n)]
             for y in range(n)]
    upper = [[0 for x in range(n)]
             for y in range(n)]
 
    # Decomposing matrix into Upper
    # and Lower triangular matrix
    for i in range(n):
 
        # Upper Triangular
        for j in range(i, n):
 
            # Summation of L(i, j) * U(j, k)
            sum = 0
            for k in range(i):
                sum += (lower[i][k] * upper[k][j])
 
            # Evaluating U(i, k)
            upper[i][j] =  M[i][j] - sum
 
        # Lower Triangular
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1  # Diagonal as 1
            else:
 
                # Summation of L(k, j) * U(j, i)
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
 
                # Evaluating L(k, i)
                lower[k][i] = ((M[k][i] - sum) /
                                  upper[i][i])
    print("Lower Triangular")

    # Displaying the result :
    for i in range(n):
        print(lower[i])

    print ("Upper Triangular")
    for i in range(n):
        print(upper[i])
    LUDoolite.L = lower
    LUDoolite.U = upper

#forward backward subsitution
def forwardbackward(L, U, b):
    swapRows(b, 0, swapnumber, len(b[0]))  #swapping rows of b accoring to partial pivoting of original matrix
    y = [[0] for r in range(len(b))] #empty list to store y
#forward substitution

    for i in range(len(b)):
        y[i][0] = b[i][0]  #storing b in empty y
        for j in range(i):
            y[i][0]=y[i][0]-(L[i][j]*y[j][0])  #formula to calculate y
   
        y[i][0] = y[i][0]/L[i][i]  
    n = len(y)

    x = [[0] for r in range(len(b))]  #empty list to store x
    if U[n-1][n-1] == 0:  #checking if diagonal elements are 0
        raise ValueError

#backward substitution
    for i in range(n-1, -1, -1):  
        x[i][0] = y[i][0]    #temporarily storing y in x
        for j in range(i+1,n):
            x[i][0] = x[i][0] -(U[i][j]*x[j][0])  #formula for x
        x[i][0] = x[i][0]/U[i][i]    

    return x
 
 #Crout algorithm to calculate Lu
def LUCrout(A):
    PartialPivot(A, 0, len(A), len(A[0])) #partial pivoting the given matrix
    n = len(A)        
    L = [[0 for x in range(n)]  #empty matrix for L and U
             for y in range(n)]
    U = [[0 for x in range(n)]
             for y in range(n)]
 
    for j in range(n):
        U[j][j] = 1             # set the j,j-th entry of U to 1
        for i in range(j, n):  # starting at L[j][j], solve j-th column of L
            a = float(A[i][j])
            for k in range(j):
                a -= L[i][k]*U[k][j]  #summating
            L[i][j] = a
        for i in range(j+1, n):# starting at U[j][j+1], solve j-th row of U
            tempU = float(A[j][i])
            for k in range(j):
                tempU -= L[j][k]*U[k][i] #summating
            U[j][i] = tempU/L[j][j]
    print ("Lower Triangular")
    for i in L:
        print (i)
    print ("\nUpper Triangular")
    for j in U:
        print (j)
    LUCrout.L = L
    LUCrout.U = U
    return L,U

#forward backward substitution to calculate inverse using LU
def forwardbackwardI(L, U, b):
    y = [[0 for c in range(len(b[0]))] for r in range(len(b))]
    for i in range(len(b)):
        for k in range (len(b[0])): #looping over the coloumns to calculate y
            y[i][k] = b[i][k]
            for j in range(i):
                y[i][k]=y[i][k]-(L[i][j]*y[j][k]) #formula for y
    
            y[i][k] = y[i][k]/L[i][i] 

    n = len(y)


    x = [[0,0,0,0] for r in range(len(b))]
    if U[n-1][n-1] == 0: #checking if diagonal elements are zero
        raise ValueError

    for i in range(n-1, -1, -1):
        for k in range (len(b[0])): #iterating over coloumns to calculate x
            x[i][k] = y[i][k]
            for j in range(i+1,n):
                x[i][k] = x[i][k] -(U[i][j]*x[j][k]) #formula for x
            x[i][k] = x[i][k]/U[i][i]    
    
    print ("The inverse of the given matrix is " )
    for i in x:
        print (i)
    return(x)

def inverseLU(M,I): #function to call inverse using LU
    if M[1][1] == 0 and M[0][1] != 0:
        swapRows(M, 0,1,4) #if diagonal element is 0, swaps to prevent 0 determinant
    LUDoolite(M)

    L = LUDoolite.L #making the given matrix into LU form
    U = LUDoolite.U 
    return forwardbackwardI(L,U,I)

#cholesky decomposition
def Cholesky_Decomposition(matrix, n):
    PartialPivot(matrix, 0, 4, 4) #partial pivoting the matrxix
 
    lower = [[0 for x in range(n)]
                for y in range(n)]
 
    # Decomposing a matrix
    # into Lower Triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0
 
            # summation for diagonals
            if (j == i):
                for k in range(j):
                    sum1 += (lower[j][k])**(2)
                lower[j][j] = round(float((matrix[j][j] - sum1)**(0.5)),2)  #formula to calculate diagonals
            else:
                 
                # Evaluating L(i, j)
                # using L(j, j)
                for k in range(j):
                    sum1 += (lower[i][k] *lower[j][k])
                if(lower[j][j] > 0):
                    lower[i][j] = round(float((matrix[i][j] - sum1) / #formula to calculate L[i][j]
                                               lower[j][j]),2)
 
    print("Lower Triangular")

    # Displaying the result :
    for i in range(n):
        print(lower[i])

    print ("Upper Triangular")
    for i in range(n):
        print(transpose(lower)[i])


    Cholesky_Decomposition.L = lower
    Cholesky_Decomposition.U = transpose(lower) #storng the U as trasnpose of L matrix
