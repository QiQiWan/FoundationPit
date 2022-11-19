import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from FoundationPit import *

# Double integral
def integrate2(f, x, y):
    g = integrate(f, x)
    g = integrate(g, y)
    return g
# Double differential
def diffn(f, x, n):
    while n != 0:
        f = diff(f, x)
        n = n - 1
    return f

def DiffIntegrate(w: Function, z, steps: int,
                  lower: int, upper: int) -> Expr:
    """Approximate integral of finite difference method
    Args:
        w (Function): The function to integrate
        z (Symbol): The variable to integrate
        steps (int): Integration steps
        upper (int): Integral upper bound
        lower (int): Integral lower bound

    Returns:
        Expr: Numerical integral expression
    """
    stepSize = (upper - lower) / steps
    sum = 0
    value = lower
    for i in range(steps):
        bottom = value
        value += stepSize
        top = value
        mid = (top + bottom) / 2
        sum += w.subs(z, value) * stepSize
    return sum
    
def P0(z, gamma, K0, D):
    """Calculate the active earth pressure 

    Args:
        z (symbol): Symbol name
        gamma (float): Bulk density of soil
        K0 (float): Active earth pressure coefficient
        D (int): Width for calculation

    Returns:
        Expr: Active earth pressure
    """
    return z * gamma * K0 * D

def Pp(z, gamma, Kp, D):
    """Calculate the passive earth pressure

    Args:
        z (symbol): Symbol name
        gamma (float): Bulk desity of soil
        Kp (float): Passive earth pressure coefficient
        D (int): Width for calculation

    Returns:
        Expr: Passive earth pressure
    """
    return z * gamma * Kp * D

def Pa(z, gamma, Ka, D):
    """Calculate the static earth pressure

    Args:
        z (symbol): Symbol name
        gamma (float): Bulk desity of soil
        Ka (float): Passive earth pressure coefficient
        D (int): With for calculation

    Returns:
        Expr: Static earth pressure
    """
    return z * gamma * Ka * D

def CalMp(z, dsa, Nsa, p0):
    """计算弯矩图

    Args:
        z (symbol): 自变量符号
        dsa (list): 支撑位置列表
        Nsa (list): 支撑支撑轴力列表
        p0 (expr): 主动土压力函数

    Returns:
        Piecewise: 弯矩分段函数 
    """
    pieces = []
    Mp = integrate2(p0, (z, 0, z), (z, 0, z))
    for i in range(len(Nsa)):
        Mp -= integrate(Nsa[i], (z, dsa[i], z))
        pieces.append((Mp, z <= dsa[i + 1]))
    Mp = Piecewise(pieces)
    return Mp

class PieceWiseFunction(object):
    """Class of the segmented function
    """
    def __init__(self, EI):
        self.periods = []
        self.funcs = []
        self.EI = EI

    def addPiece(self, func, upBoudnary):
        if len(self.periods) < 1:
            self.periods.append(upBoudnary)
            self.funcs.append(func / self.EI)
            return
        for i in range(len(self.periods)):
            if (upBoudnary < self.periods[i]):
                break
        self.periods = self.periods[:i + 1] + [upBoudnary] + self.periods[i + 1:]
        self.funcs = self.funcs[:i + 1] + [func / self.EI] + self.funcs[i + 1:]

    def getPieces(self):
        return self.periods

    def getFunctions(self):
        return self.funcs

    def getFunction(self, value):
        for i in range(len(self.periods)):
            if value <= self.periods[i]:
                return self.funcs[i]
        return False

    def subs(self, symbol, value):
        for i in range(len(self.periods)):
            if value <= self.periods[i]:
                return self.funcs[i].subs(symbol, value)
        return 0

    def CheckLimitPosi(self, symbol, lim, precisin):
        if self.subs(z, self.periods[-2]) <= lim:
            return self.periods[-2]
        posi = self.periods[-2]
        while posi <= self.periods[-1]:
            if self.subs(symbol, posi) <= lim:
                return posi
            posi += precisin
        return self.periods[-1]

def CalDeflection(z, p0, dsa, Nsa, CSet, EI):
    """Calculate the deflections of piles are calculated by given stress and stiffness of supporting pile

    Args:
        p0 (expr): Static earth pressure, a express
        dsa (list): Positions of node on the piles
        Nsa (list): List of the supporting axial force
        CSet (list): List of the coefficients of the deflection function
        EI (int): Stiffness of the piles

    Returns:
        Piecewise: Deflection function
    """
    theta0 = CSet[0]
    deflection0 = CSet[1]
    thell = []
    wl = []
    Mpr = integrate2(p0, (z, 0, z), (z, 0, z))
    index = 0
    while index < len(dsa) - 1:
        Mpr -= integrate(Nsa[index], (z, dsa[index], z))
        theta = theta0 + integrate(Mpr, (z, dsa[index], z))
        theta0 += integrate(Mpr, (z, dsa[index], dsa[index + 1]))
        deflection = -deflection0 - integrate(theta, (z, dsa[index], z))
        deflection0 += integrate(theta, (z, dsa[index], dsa[index + 1]))
        thell.append(theta)
        wl.append((-deflection / EI, z <= dsa[index + 1]))
        index += 1

    deflection = Piecewise(*wl)
    return piecewise_fold(deflection)

def CalP0ModifyFunc(z, wf, dsa, p0Lim, p0, pa):
    """Calculate the modified active earth pressure function 

    Args:
        z (symbol): Symbol
        wf (Piecewise): Deflection function
        dsa (list): Positions of node on the piles
        p0Lim (float): Active limit earth pressure displacement
        p0 (expr): Static earth pressure
        pa (expr): Active earth pressure

    Returns:
        expr: Modified active earth pressure function 
    """

    modifyP0 = piecewise_fold(p0 - (p0 - pa) * (wf / p0Lim))

    Y = [i / 10 for i in range(dsa[-1] * 10 + 1)]
    X = [float(modifyP0.subs(z, i / 10)) for i in range(dsa[-1] * 10 + 1)]
    coeffs = np.polyfit(Y, X, 5)
    P0Eq = 0
    for i in range(len(coeffs)):
        P0Eq += coeffs[i] * z ** (len(coeffs) - i - 1)

    modifyP0 = P0Eq
    return X, modifyP0

def CalPpModifyFunc(z, wf, dsa, ppLim, p0, pp):
    """Calculate the modified passive earth pressure function 

    Args:
        z (symbol): Symbol
        wf (Piecewise): Deflection function
        dsa (list): Positions of node on the piles
        ppLim (float): Passive limit earth pressure displacement
        p0 (expr): Static earth pressure
        pp (expr): Active earth pressure

    Returns:
        expr: Modified passive earth pressure function 
    """
    modifyPp = piecewise_fold(p0 + (pp - p0) * (wf / ppLim))
    Y = [i / 10 for i in range(dsa[-1] * 10 + 1)]
    X = [float(modifyPp.subs(z, i / 10)) for i in range(dsa[-1] * 10 + 1)]
    coeffs = np.polyfit(Y, X, 5)
    PpEq = 0
    for i in range(len(coeffs)):
        PpEq += coeffs[i] * z ** (len(coeffs) - i - 1)

    modifyPp = PpEq

    return X, modifyPp
#############################################
def WriteExcel(z: Symbol, wl: Expr, wr: Expr, wm: Expr, L2: int):
    dir = 'data/'
    name = f'桩长{L2}m'
    with open(f'{dir}{name}.csv', 'w+') as f:
        f.write(f'桩长{L2}m,深侧变形m,浅侧变形m\n')
        for i in range(231):
            wrl = wl.subs(z, i / 10) * 1000
            wrr = ''
            wrm = ''
            if i / 10 <= L2:
                wrr = -wr.subs(z, i / 10) * 1000
            if i / 10 <= 13: 
                wrm = -wm.subs(z, i / 10) * 1000
            f.write(f'{i / 10}, {wrl}, {wrr}, {wrm}\n')

# 定义基坑计算必须的参数
E = 200E9
h = 0.8
I = h ** 3 / 12
EI = 1.28E9
v = 0.3
phi = 37.2
c = 0
# modify the passive soil pressure
ks = (0.2 * phi ** 2 - phi + c) / 0.01
K0 = 0.5
Kp = float(tan(pi / 4 + phi / 2 / 180 * pi) ** 2)
Ka = float(tan(pi / 4 - phi / 2 / 180 * pi) ** 2)
# Ka = 0.3
gamma = 18
D = 1

supportCount = 2

L1 = L2 = 23
H1 = H2 = 12.5
L2 = 14
H2 = 6.5
B = 30

# L1, H1, L2, H2 = symbols('L_1 H_1 L_2 H_2')
Ns = [symbols(f'N_{i + 1}') for i in range(supportCount + 2)]
Nsn = {i: 0 for i in Ns}
NEA = [1.21E9 for i in Ns]
ds = [3 * i for i in range(supportCount)]
z = symbols('z', positive=True)

concreteMaterial = ConcreteMaterial('钢筋混凝土支撑', 25, E, v)
intervals = []
soilMaterial1 = SoilMaterial('粉质黏土', 19.5E3, 6.5E6, phi=27, c=49.92E3)
intervals.append({
    'top': 0,
    'bottom': 8,    
})
soilMaterial2 = SoilMaterial('中砂', 15E3, 25E6, phi=31.6, c=0)
intervals.append({
    'top': 8,
    'bottom': 17,    
})
soilMaterial3 = SoilMaterial('砾砂', 20E3, 32E6, phi=35, c=0)
intervals.append({
    'top': 17,
    'bottom': 23,    
})
boreHole = BoreHole([soilMaterial1, soilMaterial2, soilMaterial3], intervals)
supportMaterial = SupportMaterial('钢混支撑', 25, E, 0.000927)

leftWall = UndergroundDiaphragmWall(L1, H1, concreteMaterial)
rightWall = UndergroundDiaphragmWall(L2, H2, concreteMaterial)

supports = [HorizontalSupport(supportMaterial, 30) for _ in range(supportCount)]

foundationPit = FoundationPit(leftWall, rightWall, H1, H2, supports,
                              supportCount, ds, B, D, boreHole)

#####################################################

# 构造位移函数，形函数为多项式
dsal = [0] + ds + [H1] + [L1]
dsar = [0] + ds + [H2] + [L2]
dsam = [0] + ds + [6.5] + [13]

def EnergyFormulate(z, foundationPit: FoundationPit, L2=16, hm=8):
    # 计算坑中坑土压力参数
    PalimWm = 0.005 * 13
    PplimWm = 0.05 * 13
    PalimWm

    # def MinPotentialEnergy(foundationPit: FoundationPit, hm=6, L2=16):
    wl = 0 # 位移函数，是一个一元高次多项式
    wr = 0 # 位移函数，是一个一元高次多项式
    wm = 0 # 坑中坑位移函数
    hm = 8
    dsar[-1] = L2

    CSet = [] # wl多项式参数列表
    BSet = [] # wr多项式参数列表
    DSet = [] 

    index = 1
    # 循环里面创建位移多项式
    for j in range(hm):
        Cp = symbols('C_%d' % index)
        wl += Cp * z ** j
        Bp = symbols('B_%d' % index)
        wr += Bp * z ** j
        Dp = symbols('D_%d' % index)
        wm += Dp * z ** j
        index += 1
        CSet.append(Cp)
        BSet.append(Bp)
        DSet.append(Dp)

    # 下面的等式是的位移函数满足底端变形为0的边界条件
    wl = wl * (z - dsal[-1]) ** 2
    wr = wr * (z - dsar[-1]) ** 2
    wm = wm * (z - dsam[-1]) ** 2

    # 将两个位移函数用支撑变形协调条件连接起来
    eqs = []
    for i in range(supportCount):
        eqs.append(Eq(wl.subs(z, ds[i]) + wr.subs(z, ds[i]), Ns[i] * B / NEA[i]))
        eqs.append(Eq(wl.subs(z, 6.5 + ds[i]) + wm.subs(z, ds[i]), Ns[i + 2] * B / 2 / NEA[i + 2]))

    result = solve(eqs, BSet + DSet)
    wr = wr.subs(result)
    wm = wm.subs(result)
    for key in result:
        if key in BSet:
            BSet.remove(key)
        if key in DSet:
            DSet.remove(key)

    # 支护桩外力做功
    Wl = 0 # 左侧支护桩的外力功
    Wr = 0 # 右侧支护桩的外力功
    Wm = 0 # 坑中坑支撑外力做功

    # steps = 100
    # 主动区外力做功
    # sl, sr = symbols('s_l s_r')
    # WlTemp = integrate(foundationPit.Pa(z, sl, LRSide.LeftSide), (sl, 0, wl))
    # Wl += DiffIntegrate(WlTemp, z, steps, 0, dsal[-1])
    # WrTemp = integrate(foundationPit.Pa(z, sr, LRSide.RightSide), (sr, 0, wr))
    # Wr += DiffIntegrate(WrTemp, z, steps, 0, dsar[-1])
    Pacrl = foundationPit.Pacr(z, LRSide.LeftSide)
    P0l = foundationPit.P0(z, LRSide.LeftSide)
    Pacrr = foundationPit.Pacr(z, LRSide.RightSide)
    P0r = foundationPit.P0(z, LRSide.RightSide)
    Wl += integrate((2 * P0l + (Pacrl - P0l) * 1.9 * wl / foundationPit.PalimWl) * wl/ 2, (z, 0, dsal[-1]))
    Wr += integrate((2 * P0l + (Pacrr - P0r) * 1.9 * wr / foundationPit.PalimWr) * wr/ 2, (z, 0, dsar[-1]))
    Wm += integrate((2 * P0l + (Pacrl - P0l) * 1.9 * wm / PalimWm) * wm / 2, (z, 0, dsam[-1]))

    # 被动区外力做功
    Ppcrl = foundationPit.Ppcr(z - H1, LRSide.LeftSide)
    P0l = foundationPit.P0(z - H1, LRSide.LeftSide)
    Ppcrr = foundationPit.Ppcr(z - H2, LRSide.RightSide)
    P0r = foundationPit.P0(z - H2, LRSide.RightSide)
    Ppcrm = foundationPit.Ppcr(z - dsam[-2], LRSide.LeftSide)
    P0m = foundationPit.P0(z - dsam[-2], LRSide.LeftSide)
    # WlTemp = integrate(foundationPit.Pp(z, sl, LRSide.LeftSide), (sl, 0, wl))
    # Wl -= DiffIntegrate(WlTemp, z, steps, dsal[-2], dsal[-1])
    # WrTemp = integrate(foundationPit.Pp(z, sr, LRSide.RightSide), (sr, 0, wr))
    # Wr -= DiffIntegrate(WrTemp, z, steps, dsar[-2], dsar[-1])
    Wl -= integrate((2 * P0l + (Ppcrl - P0l) * 1.9 * wl / foundationPit.PplimWl) * wl / 2, (z, dsal[-2], dsal[-1]))
    Wr -= integrate((2 * P0r + (Ppcrr - P0r) * 1.9 * wr / foundationPit.PplimWr) * wr / 2, (z, dsar[-2], dsar[-1]))
    Wm -= integrate((2 * P0m + (Ppcrm - P0m) * 1.9 * wm / PplimWm) * wm / 2, (z, dsam[-2], dsam[-1]))

    # Wl -= integrate(((2 * P0(z - dsal[-2], gamma, K0, D) + (Pp(z - dsal[-2], gamma, Kp, D)\
    # - P0(z - dsal[-2], gamma, K0, D)) * wl / PplimWl) * wl / 2), (z, dsal[-2], dsal[-1]))
    # Wr -= integrate(((2 * P0(z - dsar[-2], gamma, K0, D) + (Pp(z - dsar[-2], gamma, Kp, D)\
    # - P0(z - dsar[-2], gamma, K0, D)) * wr / PplimWr) * wr / 2), (z, dsar[-2], dsar[-1]))

    # 支护桩变形能
    Ul = 0 # 左侧支护桩的应变能
    Ur = 0 # 右侧支护桩的应变能
    Um = 0 # 坑中坑支护桩应变能
    EIl = foundationPit.LeftWall.EI
    EIr = foundationPit.RightWall.EI
    Ul += integrate(EI * wl.diff(z, 2) ** 2 / 2, (z, 0, dsal[-1]))
    Ur += integrate(EI * wr.diff(z, 2) ** 2 / 2, (z, 0, dsar[-1]))
    Um += integrate(EI * wm.diff(z, 2) ** 2 / 2, (z, 0, dsam[-1]))


    # # 铁木辛科梁剪切应变能
    # alpha = 6 / 5
    # GAl = foundationPit.LeftWall.G * foundationPit.LeftWall.h * foundationPit.D
    # GAr = foundationPit.RightWall.G * foundationPit.LeftWall.h * foundationPit.D
    # Ul += alpha * EI ** 2 / 2 / GAl * integrate(wl.diff(z, 3) ** 2, (z, 0, dsal[-1]))
    # Ur += alpha * EI ** 2 / 2 / GAr * integrate(wr.diff(z, 3) ** 2, (z, 0, dsar[-1]))
    # Um += alpha * EI ** 2 / 2 / GAr * integrate(wm.diff(z, 3) ** 2, (z, 0, dsam[-1]))


    # 支撑轴力应变能
    Ne = 0
    for i in range(supportCount):
        Ne += Ns[i] ** 2 * B / 2 / NEA[i]
        Ne += Ns[i + 2] ** 2 * B / 4 / NEA[i]
    TotU = Ul + Ur + Um - Wl - Wr - Wm + Ne
    eqs = []
    # 最小势能原理，计算总势能变分
    for i in BSet:
        eqs.append(Eq(TotU.diff(i), 0))
    for i in CSet:
        eqs.append(Eq(TotU.diff(i), 0))
    for i in DSet:
        eqs.append(Eq(TotU.diff(i), 0))
    for i in Ns:
        eqs.append(Eq(TotU.diff(i), 0))
    result = solve(eqs, DSet + BSet + CSet + Ns)
    
    return result, wl, wr, wm

result, wl, wr, wm = EnergyFormulate(z, foundationPit=foundationPit, L2=16)

# 绘图查看计算结果
# head = dsal[-2] * 10
head = 0
count = L1 * 10 + 1
# count = int(foundationPit.PplimWl * 1000 + 1)

w = symbols('w')
plt.ylim(int(count / 10), head / 10)
# func = (2 * P0(z - dsar[-2], gamma, K0, D) + (Pp(z - dsar[-2], gamma, Kp, D)\
#     - P0(z - dsar[-2], gamma, K0, D)) * wr.subs(result) / PplimWr)
# func = (2 * p0 - (p0 - pa) * wr.subs(result) / P0limWr)
# func = wr.subs(result)

# func = ((Pp(z - dsar[-2], gamma, Kp, D)\
#     - P0(z - dsar[-2], gamma, K0, D)) * wl.subs(result) / PplimWr)
# func = P0(z - dsal[-2], gamma, K0, D) + ((Pp(z - dsal[-2], gamma, Kp, D)\
#     - P0(z - dsal[-2], gamma, K0, D)) * wl.subs(result) / PplimWr)
deepth = 16

# func = foundationPit.Pp(z, w, LRSide.LeftSide).subs(z, deepth)

func = wl.subs(result).diff(z, 2) * EI / 1E3
Pacrl = foundationPit.Ppcr(z - H1, LRSide.LeftSide).subs(z, deepth)
P0l = foundationPit.P0(z - H1, LRSide.LeftSide).subs(z, deepth)

# with open('result.csv', 'w+') as f:
#     f.write('deformation, origin, fitness\n')
#     for i in range(head, count):
#         f.write(f'{i / 1000}, {func.subs(w, i / 1000)}, {(P0l + (Pacrl - P0l) * 1.9 * w / foundationPit.PplimWl).subs(w, i / 1000)}\n')

X = [func.subs(z, i / 10) for i in range(head, count)]
print(max(X))

plt.plot(X, [i / 10 for i in range(head, count)])
# plt.plot([i / 1000 for i in range(head, count)], [func.subs(w, i / 1000) for i in range(head, count)])
# plt.plot([i / 1000 for i in range(head, count)], [(P0l + (Pacrl - P0l) * 1.9 * w / foundationPit.PplimWl).subs(w, i / 1000) for i in range(head, count)])

plt.savefig('fig.png', dpi=1000)
plt.show()