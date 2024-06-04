import numpy as np
import math as ma
import cmath as cma
import sympy as sym
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


class compound(object):#化合物类，用于读取数据表中的数据
    def __init__(self,s,metal=False):#如果材料是金属，则会拥有消光系数，记得标上metal=True
        self.name = s
        self.compoundfile = open(self.name+"-Refractive Index.csv","r",encoding="utf-8").readlines()#需要数据表文件在同一个目录下面
        self.metal=metal
        j = 0
        self.wavelength_list = []
        self.reflective_list = []
        self.extinction_list = []
        for i in self.compoundfile:
            self.compoundfile[j] = i.replace("\n", "")
            j += 1
        for i in range(len(self.compoundfile)):
            if i == 0:
                continue
            k = self.compoundfile[i].split(",")
            self.wavelength_list.append(float(k[0]))
            self.reflective_list.append(float(k[1]))
            if self.metal:#读取金属的消光系数
                self.extinction_list.append(float(k[2]))
        if self.metal:
            for i in range(len(self.reflective_list)):
                self.reflective_list[i] = self.reflective_list[i]+self.extinction_list[i]*1j#计算并存入金属的复折射率

class grid(object):#网格类，用于读取网格划分信息
    def __init__(self):#用点的字典绑定点的序数和坐标
        self.point_index = open("point file","r",encoding="utf-8").readlines()
        for i in self.point_index:
            self.point_index[j] = i.replace("\n", "")
            j += 1
class material(object):#TMM矩阵计算
    def __init__(self, layerlist, thicklist,refraction_index,theta_air,k):#输入材料列表与厚度和入射角度
        self.layerlist = layerlist
        self.thicklist = thicklist        
        self.refraction_index = refraction_index #需要按顺序写下每一层指定波长下的折射率
        self.k = k  #n为复数时，k也是复数
        self.theta_air = theta_air
    
      
    def calcu_thetalist(self,theta_air):#设置为字典可以存储复数
        self.thetalist = {}
        for i in range(len(self.refraction_index)):
            if i == 0 :
                self.thetalist[i] = cma.asin(cma.sin(theta_air)/self.refraction_index[i])
            else:
                self.thetalist[i] = cma.asin(cma.sin(self.thetalist[i-1])*self.refraction_index[i-1]/self.refraction_index[i])
        return self.thetalist
    
    def cal_beta(self,theta_air,k):#计算传播方向光线波速的分量，每层都不相同,k与w的关系是k=n*w/c,空气中n为1，这是真空中的
        return k*ma.sin(theta_air)
    
    def cal_kalist(self,k,theta_air):#计算垂直传播方向光线波速的分量，每层都不同，在角度引入复数
        self.kalist = {}
        thetalist = self.calcu_thetalist(theta_air)
        for i in range(len(self.refraction_index)):
            self.kalist[i] = k*self.refraction_index[i]*cma.cos(thetalist[i])
        return self.kalist
    def Da_smatrix(self,theta_air):#计算每层的D矩阵和它的逆矩阵，放在同一个层数中，s波
        thetalist = self.calcu_thetalist(theta_air)
        shape = (len(self.refraction_index),2,2,2)
        Da_smatrix = np.zeros(shape, dtype=complex)
        for i in range(len(self.refraction_index)):
            Da_smatrix[i][0] = np.array([[1,1],[self.refraction_index[i]*cma.cos(thetalist[i]),\
                                             -1*self.refraction_index[i]*cma.cos(thetalist[i])]])
            Da_smatrix[i][1] = np.linalg.inv(Da_smatrix[i][0])
        return Da_smatrix
    def Da_pmatrix(self,theta_air):#计算每层的D矩阵和逆矩阵，p波
        thetalist = self.calcu_thetalist(theta_air)
        shape = (len(self.refraction_index),2,2,2)
        Da_pmatrix = np.zeros(shape, dtype=complex)
        for i in range(len(self.refraction_index)):
            Da_pmatrix[i][0]= np.array([[cma.cos(self.thetalist[i]),cma.cos(thetalist[i])],\
                                      [self.refraction_index[i],-1*refraction_index[i]]])
            Da_pmatrix[i][1]= np.linalg.inv(Da_pmatrix[i][0])
        return Da_pmatrix
    def Pa_matrix(self,k,theta_air):#计算每层的传播矩阵,逆矩阵不存在
        kalist = self.cal_kalist(k,theta_air)
        shape = (len(self.refraction_index),2,2)
        Pa_matrix = np.zeros(shape, dtype=complex)
        shape1 = (len(self.refraction_index),1)
        fhia = np.zeros(shape1, dtype=complex)
        for i in range(len(self.refraction_index)):
            fhia[i] = kalist[i]*self.thicklist[i]
            Pa_matrix[i] = np.array([[cma.exp(1j*fhia[i]),0],\
                                      [0,cma.exp(-1j*fhia[i])]])
        return Pa_matrix
    def cal_airmatrix(self,theta_air,distance,k):#这里考虑从发光点到材料的入射点，distance是两者的垂直距离

        D_airmatrixs = np.array([[1,1],[cma.cos(theta_air),\
                                             -cma.cos(theta_air)]])
        D_airmatrixp = np.array([[cma.cos(theta_air),cma.cos(theta_air)],\
                                      [1,-1]])
        Pa_airmatrix = np.array([[cma.exp(1j*distance*k*cma.sin(theta_air)),0],\
                                      [0,cma.exp(-1j*distance*k*cma.sin(theta_air))]])
        return D_airmatrixs,D_airmatrixp,Pa_airmatrix
    def M_matrix_s(self,theta_air,distance,k):#总体的传播矩阵，s波，从入射材料的点算起，到最后一层,n+1层的真空处出来
        airmatrix = self.cal_airmatrix(theta_air,distance,k)
        Da_smatrix = self.Da_smatrix(theta_air)
        Pa_matrix = self.Pa_matrix(k,theta_air)
        if len(self.layerlist) == 1:
            T = np.linalg.inv(airmatrix[0])@Da_smatrix[0][0]@Pa_matrix[0]@Da_smatrix[0][1]@airmatrix[0]
        else:
            for i in range(len(self.layerlist)):
                if i == 0 :
                    T = np.linalg.inv(airmatrix[0])@Da_smatrix[0][0]@Pa_matrix[0]@Da_smatrix[0][1]@Da_smatrix[1][0]
                elif i < len(self.layerlist-1):
                    T = T@Da_smatrix[i-1][1]@Da_smatrix[i][0]@Pa_matrix[i]@Da_smatrix[i][1]@Da_smatrix[i+1][0]
                else:
                    T= T@Da_smatrix[i-1][1]@Da_smatrix[i][0]@Pa_matrix[i]@Da_smatrix[i][1]@airmatrix[0]
        
        return T
    def M_matrix_p(self,theta_air,distance,k):#总体的传播矩阵，p波，从入射材料的点算起，到最后一层,n+1层的真空处出来
        airmatrix = self.cal_airmatrix(theta_air,distance,k)
        Da_pmatrix = self.Da_smatrix(theta_air)
        Pa_matrix = self.Pa_matrix(k,theta_air)
        if len(self.layerlist) == 1:
            T = np.linalg.inv(airmatrix[0])@Da_pmatrix[0][0]@Pa_matrix[0]@Da_pmatrix[0][1]@airmatrix[0]
        else:
            for i in range(len(self.layerlist)):
                if i == 0 :
                    T = np.linalg.inv(airmatrix[0])@Da_pmatrix[0][0]@Pa_matrix[0]@Da_pmatrix[0][1]@Da_pmatrix[1][0]
                elif i < len(self.layerlist-1):
                    T = T@Da_pmatrix[i-1][1]@Da_pmatrix[i][0]@Pa_matrix[i]@Da_pmatrix[i][1]@Da_pmatrix[i+1][0]
                else:
                    T= T@Da_pmatrix[i-1][1]@Da_pmatrix[i][0]@Pa_matrix[i]@Da_pmatrix[i][1]@airmatrix[0]
        
        return T
    def M_every_matrix_s(self,theta_air,distance,k,n):#第n层(从1开始)从进入到右层边界处同种材料内的位置的矩阵,s波
        airmatrix = self.cal_airmatrix(theta_air,distance,k)
        Da_smatrix = self.Da_smatrix(theta_air)
        Pa_matrix = self.Pa_matrix(k,theta_air)
        shape = (2,2)
        Ts = np.zeros(shape, dtype=complex)
        for i in range(n):
            if i == 0 :
                Ts = np.linalg.inv(airmatrix[0])@Da_smatrix[0][0]@Pa_matrix[0]
            elif i < n or i == n:
                Ts = Ts@Da_smatrix[i-1][1]@Da_smatrix[i][0]@Pa_matrix[i]
            
        return Ts
    def M_every_matrix_p(self,theta_air,distance,k,n):#第n层从进入到右层边界处同种材料内的位置的矩阵,p波
        airmatrix = self.cal_airmatrix(theta_air,distance,k)
        Da_pmatrix = self.Da_pmatrix(theta_air)
        Pa_matrix = self.Pa_matrix(k,theta_air)
        shape = (2,2)
        Tp = np.zeros(shape, dtype=complex)
        for i in range(n):
            if i == 0 :
                Tp = np.linalg.inv(airmatrix[0])@Da_pmatrix[0][0]@Pa_matrix[0]
            elif i < n or i == n:
                Tp = Tp@Da_pmatrix[i-1][1]@Da_pmatrix[i][0]@Pa_matrix[i]
            
        return Tp
    


    def epsilon(self,n,wavelength):#计算具有w波长的光通过第n层材料的光程,垂直情况下
        return 2*(ma.pi/wavelength)*self.refraction_index(n,wavelength)*self.thicklist[n]/1000
    def cos_epsilon(self,n,w):#复三角函数计算，防止数据溢出
        return (cma.cos(self.epsilon(n,w).real))*(cma.cosh(self.epsilon(n,w).imag))-(cma.sin(self.epsilon(n,w).real))*(cma.sinh(self.epsilon(n,w).imag))*1j
    def sin_epsilon(self,n,w):#复三角函数计算，防止数据溢出
        return (cma.sin(self.epsilon(n, w).real)) * (cma.cosh(self.epsilon(n, w).imag)) + (cma.cos(self.epsilon(n, w).real)) * (cma.sinh(self.epsilon(n, w).imag)) * 1j
    def eigenmatrix_s(self,n,w):#单层s波特征矩阵,具体计算方法见论文
        return np.array([[self.cos_epsilon(n,w),(-1j/self.refraction_index(n,w))*self.sin_epsilon(n,w)],[(-1j*self.refraction_index(n,w))*self.sin_epsilon(n,w),self.cos_epsilon(n,w)]])
    def eigenmatrix_p(self,n,w):#单层p波特征矩阵
        return np.array([[self.cos_epsilon(n,w),(-1j*self.refraction_index(n,w))*self.sin_epsilon(n,w)],[(-1j/self.refraction_index(n,w))*self.sin_epsilon(n,w),self.cos_epsilon(n,w)]])
    def T_matrix_s(self,w):#s波特征矩阵
        for i in range(len(self.layerlist)):
            if i== 0:
                T = self.eigenmatrix_s(i,w)
            else:
                T =np.dot(T,self.eigenmatrix_s(i,w))
        return T
    def T_matrix_p(self,w):#p波特征矩阵
        for i in range(len(self.layerlist)):
            if i== 0:
                T = self.eigenmatrix_p(i,w)
            else:
                T =np.dot(T,self.eigenmatrix_p(i,w))
        return T
    def t(self,w):#s波透射系数
        T = self.T_matrix_s(w)
        return 2/(T[0,0]+T[0,1]+T[1,1]+T[1,0])
    def r(self,w):#s波反射系数
        T = self.T_matrix_s(w)
        return (T[0,0]+T[0,1]-T[1,0]-T[1,1])/(T[0,0]+T[0,1]+T[1,1]+T[1,0])
    def t_p(self,w):#p波透射系数
        T = self.T_matrix_p(w)
        return 2/(T[0,0]+T[0,1]+T[1,1]+T[1,0])
    def r_p(self,w):#p波反射系数
        T = self.T_matrix_p(w)
        return (T[0,0]+T[0,1]-T[1,0]-T[1,1])/(T[0,0]+T[0,1]+T[1,1]+T[1,0])
    def R(self,w):#发射比计算
        return (abs(self.r(w))**2+abs(self.r_p(w))**2)/2
    def T(self,w):#透射比计算
        return (abs(self.t(w))**2+abs(self.t_p(w))**2)/2
    def radiation(self,w):#发射光谱系数的计算
        return 1-self.R(w)-self.T(w)

def cal_every_edge(mater,A,B,theta_air,distance,k,thicklist,type):
    shape = (len(thicklist),2)
    list2 = np.zeros(shape, dtype=complex)
    if type == 'Es':
        for i in range(len(thicklist)):
            T = mater.M_every_matrix_s(theta_air,distance,k,i+1)
            list2[i]= np.array(cal_field2(A,B,T))
    elif type == 'Ep':
        for i in range(len(thicklist)):
            T = mater.M_every_matrix_p(theta_air,distance,k,i+1)
            list2[i]= cal_field2(A,B,T)
    
    return list2

def cal_point(x,thicklist,kalist,list1,refraction_index):
    #判断在第几层，从1开始
    thickness = []

    for i in range(len(thicklist)):
        if i ==0:
            thickness.append(thicklist[i]) 
        else:
            thickness[i] = thickness[i-1]+thicklist[i]
    for i in range(len(thickness)): 
        if x <= thickness[i]:
            xj = i
            break
        else:
            i += 1
    Ai = list1[xj][0]
    Bi = list1[xj][1]
    dx =  thickness[xj] - x
    kax = kalist[xj]
    n = refraction_index[xj]
    T = np.array([[np.exp(1j*dx*kax),0],[0,np.exp(-1j*dx*kax)]])
    Ax,Bx = cal_field3(Ai,Bi,T)
    strenth = intensity(Ax,Bx,n)
    return strenth

def cal_field(A,M):#返回的是一个元组，A是输入光线的场强，B是射出，最终射出向左的场强设置为0，解出的An是最终射出的向右的场强，B是初始的向左的场强
        An1 = sym.symbols('An1')
        B = sym.symbols('B')
        Eq1 = M[0][0] * An1-A
        Eq2 = M[1][0] * An1-B 
        dict1 = sym.solve([Eq1,Eq2],[An1,B])
        An1 = dict1[An1]
        B = dict1[B] 
        An1 = complex(An1)
        B =complex(B)
        return An1,B

def cal_field2(A,B,T):
    Ai = sym.symbols('Ai')
    Bi = sym.symbols('Bi')
    Eq1 = T[0][0]* Ai + T[0][1]* Bi - A
    Eq2 = T[1][0]* Ai + T[1][1]* Bi - B
    dict1 = sym.solve([Eq1,Eq2],[Ai,Bi])
    Ai = dict1[Ai]
    Bi = dict1[Bi] 
    Ai = complex(Ai)
    Bi =complex(Bi)
    return Ai,Bi

def cal_field3(Ai,Bi,T):
    Ax = sym.symbols('Ax')
    Bx = sym.symbols('Bx')
    Eq1 = T[0][0]* Ax + T[0][1]* Bx - Ai
    Eq2 = T[1][0]* Ax + T[1][1]* Bx - Bi
    dict1 = sym.solve([Eq1,Eq2],[Ax,Bx])
    Ax = dict1[Ax]
    Bx = dict1[Bx]
    Ax = complex(Ax)
    Bx =complex(Bx)
    return Ax,Bx

def cal_n(A,B,C,wavelength):#利用柯西公式计算折射率
    return A+B/wavelength**2+C/wavelength**3

def intensity(A,B,n):
    u0 = 4*np.pi*1e-7
    e0 = 8.854*1e-12
    C = (A+B).real**2+(A+B).imag**2
    return C*(1/2)*np.sqrt(e0/u0)*n.real

def color_map(data, cmap):
    """数值映射为颜色"""
    
    dmin, dmax = np.nanmin(data), np.nanmax(data)
    cmo = plt.cm.get_cmap(cmap)
    cs, k = list(), 256/cmo.N
    
    for i in range(cmo.N):
        c = cmo(i)
        for i in range(int(i*k), int((i+1)*k)):
            cs.append(c)
    cs = np.array(cs)
    data = np.uint8(255*(data-dmin)/(dmax-dmin))
    
    return cs[data]

if __name__ == "__main__":
    #需要场强的分布时，不考虑损耗，只需要计算出每层的A,B
    #这里需要统一物理单位，距离统一为微米，角度和折射率无单位，波速k的单位是1/微米,电场强度与A设置的相同
    #给出初始As,Bs的值和Ap,Bp的值分别计算，材料本身性质已经确定。
    Degrees = 0
    theta_air = Degrees*ma.pi/180 #角度转化为弧度
    distance= 10
    layerlist = ['Si']
    thicklist = [10]
    refraction_index = [1.2]
    k = 0.2
    mater = material(layerlist,thicklist,refraction_index,theta_air,k)
    M = mater.M_matrix_s(theta_air,distance,k)
    #print(mater.calcu_thetalist(theta_air))
    #print(mater.cal_kalist(k,theta_air))
    #print(mater.Da_smatrix(theta_air))
    #print(mater.Pa_matrix(k,theta_air))
    #print(M)
    kalist = mater.cal_kalist(k,theta_air)
    A = 3 
    An1,B = cal_field(A,M)
    list1 = cal_every_edge(mater,A,B,theta_air,distance,k,thicklist,type='Es')
    x = np.linspace(0,10,100)

    y = np.zeros(len(x), dtype=float)
    for i in range(len(x)):
        y[i] = cal_point(x[i],thicklist,kalist,list1,refraction_index)


    #x = np.linspace(0, 2*np.pi, 10)
    #y = np.sin(x)
    x1 = np.zeros(100)
    ps = np.stack((x,x1), axis=1)
    segments = np.stack((ps[:-1], ps[1:]), axis=1)

    cmap = 'viridis' # jet, hsv等也是常用的颜色映射方案
    #colors = color_map(np.cos(x)[:-1], cmap)
    colors = color_map(y[:-1], cmap)
    line_segments = LineCollection(segments, colors=colors, linewidths=200, linestyles='solid', cmap=cmap)

    fig, ax = plt.subplots()
    ax.set_xlim(np.min(x)-0.1, np.max(x)+0.1)
    ax.set_ylim(np.min(y)-0.1, np.max(y)+0.1)
    ax.add_collection(line_segments)
    cb = fig.colorbar(line_segments, cmap='jet')

    plt.show()