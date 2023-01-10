from configs import *


def tfo(z,theta):
    return 2/(z*(1-z)*e*theta**2)


def sigma(x,y,t):
    sig=qhat/(2*CF)*np.dot(x(t)-y(t),x(t)-y(t))
    return sig



def sigmalNc(x,y,t):
    sig=qhat/(Nc)*np.dot(x(t)-y(t),x(t)-y(t))
    return sig


def comb(full,last):
    current=[]
    n=len(full[0])
    for perm in last:
        for i in range(n):
            for j in range(n):
                a=[i,j]
                aset=set(a)
                if len(a)==len(aset):
                    lol=perm.copy()
                    lol[i],lol[j]=lol[j],lol[i]
                    if lol not in full:
                        current.append(lol)
                        full.append(lol)
    return current



# diaonal entries
def dia(t,z,Z):
    # The C is the constant one that is common for all the diaonals
    C=0
    for k in range(w):
        for l in range(w):
            if l > k:
                C=C+sigma(z[k],z[l],t)+sigma(Z[k],Z[l],t)-sigma(z[k],Z[l],t)-sigma(Z[k],z[l],t)
    for k in range(w):
        C=C-sigma(z[k],Z[k],t)
    
    # Cvar is not the same for all terms
    i=0
    for perm in perms:
        Cvar=0
        for k in range(w):
            Cvar=Cvar+sigma(Z[k],perm[k],t)
        # This gives the diaonal terms
        mat[i,i]= -1/4*(Nc*Cvar+1/Nc*C)
        i=i+1




# diaonal entries
def diaNc(t,z,Z):
    # Cvar is not the same for all terms
    i=0
    for perm in perms:
        Cvar=0
        for k in range(w):
            Cvar=Cvar+sigmalNc(Z[k],perm[k],t)
        # This gives the diaonal terms
        mat[i,i]= -1/4*(Nc*Cvar)
        i=i+1




# Non-diaonal entries
def nondia(t,Z): 
    k=0
    for perm in perms:
        for i in range(w):
            for j in range(w):
                if i!=j and i<j:
                    lol=perm.copy()
                    lol[i],lol[j]=lol[j],lol[i]
                    mat[k,nums[str(lol)]]=-1/4*(sigma(Z[i],perm[j],t)+sigma(Z[j],perm[i],t)-sigma(perm[i],perm[j],t)-sigma(Z[i],Z[j],t))
        k=k+1



# Non-diaonal entries
def nondiaNc(t,Z):
    k=0
    for perm in perms:
        for i in range(w):
            for j in range(w):
                if i!=j and i<j:
                    lol=perm.copy()
                    lol[i],lol[j]=lol[j],lol[i]
                    mat[k,nums[str(lol)]]=-1/4*(sigmalNc(Z[i],perm[j],t)+sigmalNc(Z[j],perm[i],t)-sigmalNc(perm[i],perm[j],t)-sigmalNc(Z[i],Z[j],t))
        k=k+1
    
    # Try to set all the entries to the right of the diagonal to zero
    for row in range(wfac):
        for col in range(wfac):
            if col > row:
                mat[row,col]=0




# Does the matrix multiplication and returns the differential equation in the large-Nc
def model(u,t):
    lis=[]
    for i in range(wfac):
        lis.append(u[i])
    v=np.array(lis)
    
    dia(t,r,R)
    
    nondia(t,R)
    #print(mat)
        
    diff=np.matmul(mat,v)
    return diff



# Does the matrix multiplication and returns the differential equation in the large-Nc
def modelNc(u,t):
    lis=[]
    for i in range(wfac):
        lis.append(u[i])
    v=np.array(lis)
    
    diaNc(t,r,R)
    
    nondiaNc(t,R)
    #print(mat)
        
    diff=np.matmul(mat,v)
    return diff



# Does the matrix multiplication and returns the differential equation

def modeldia(u,t):
    lis=[]
    for i in range(wfac):
        lis.append(u[i])
    v=np.array(lis)
    
    diaNc(t,r,R)
        
    diff=np.matmul(mat,v)
    return diff


def main(z,Z):
    global w
    w=len(z)
    # Make a matrix full of zeros that we will fill out
    global wfac
    wfac=mt.factorial(w)
    size=(wfac,wfac)
    global mat
    mat=np.zeros(size)
    
    # Make all the K! possible permutations  of the 
    full=[z]
    last=[z]
    count=[]
    while len(last)>0:
        count.append(len(last))
        last=comb(full,last)
    global perms
    perms=np.array(full)
    
    # Numbering the K! different permutations
    global nums
    nums={}
    i=0
    for perm in perms:
        nums[str(perm)]=i
        i=i+1
        
    #initial condition
    u0=[]
    j=0
    for i in count:
        u0=u0+i*[Nc**(w-j)]
        j=j+1

    # time points
    t = np.linspace(t2,t_max)

    # solve ODE in large-Nc and exact
    u = odeint(model,u0,t)
    uNc = odeint(modelNc,u0,t)
    
    mat=np.zeros(size)
    udia = odeint(modeldia,u0,t)


    np.savetxt(f'{filename}_t.csv', u, delimiter=',')
    np.savetxt(f'{filename}_full.csv', u, delimiter=',')
    np.savetxt(f'{filename}_Nc.csv', uNc, delimiter=',')
    np.savetxt(f'{filename}_dia.csv', udia, delimiter=',')