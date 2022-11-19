from enum import Enum
from typing import List

class HorizontalSupport(object):
    """Class of the horizontal support of the foundation pit
    """
    def __init__(self, material=None, spaceLength=None) -> None:
        self.Material = material
        self.SpaceLength = spaceLength
        # self.N: Axis force of the support
        self.N = None
    
    def SetMeterial(self, material):
        """Modefy the mateiral of the support

        Args:
            meterial (Material): Target material object
        """
        self.Material = material

    def SetSpaceLength(self, spaceLength):
        """Modify the calculated length
        """
        self.SpaceLength = spaceLength

class Material(object):
    """The base class of materials"""
    def __init__(self, name, gamma, E) -> None:
        self.name = name
        self.gamma = gamma
        self.E = E

class ConcreteMaterial(Material):
    """Class of concrete material, inheriate from materials class

    Args:
        Material (Material): Base class
    """
    def __init__(self, name, gamma, E, v, G = None) -> None:
        super().__init__(name, gamma, E)
        self.name = name
        self.gamma = gamma
        self.E = E
        self.v = v
        if G:
            self.G = G
        else:
            self.G = E / (2 * (1 + v))

class SupportMaterial(Material):
    """Class of support material, inheriate from materials class"""
    def __init__(self, name, gamma, E, A) -> None:
        super().__init__(name, gamma, E)
        self.name = name
        self.gamma = gamma
        self.E = E
        self.A = A
        self.EA = E * A

class EarthType(Enum):
    Rankine = 1
    Coulomb = 2

class SoilMaterial(Material):
    """Class of soil material, inheriate from materials class

    Args:
        Material (Material): Base class
    """
    def __init__(self, name: str, gamma: float,
                 E: int, phi: float, c = 0, sandy=True,
                 varepsilon=0, varphi=0, alpha=0) -> None:
        """Define the soil parameter

        Args:
            name (str): Name of the soil
            gamma (float): Weight of the soil, unit: N/m^-3
            E (int): Elastic Modulus, unit: N/m^2
            phi (float): Internal friction angle, unit: °\degree
            c (int, optional): Cohesive force, unit: Pa. Defaults to 0.
            sandy (bool, optional): Is sandy soil? Defaults to True.
            varepsilon (int, optional): Angle between the back of the retaining wall and the vertical line. Defaults to 0.
            varphi (int, optional): Friction angle of wall and soil. Defaults to 0.
            alpha (int, optional): Angle between soil layer and horizontal plane. Defaults to 0.
        """
        from math import pi, sin
        super().__init__(name, gamma, E)
        self.name = name
        self.gamma = gamma
        self.E = E
        self.phi = phi
        self.phid = phi
        self.c = c
        self.varepsilon = varepsilon
        self.varphi = varphi
        self.alpha = alpha
        self.K0 = 0
        if sandy:
            self.K0 = 1 - sin(phi / 180 * pi)
        else:
            self.K0 = 0.95 - sin(phi / 180 * pi)
        self.Ka = self.GetKa(phi, EarthType.Rankine)
        self.Kp = self.GetKp(phi, EarthType.Rankine)
        
    def GetKa(self, phi, type: EarthType):
        from math import sin, cos, tan, pi, sqrt
        phi = phi * pi / 180
        if type == EarthType.Rankine:
            return tan(pi / 4 - phi / 2) ** 2
        if type == EarthType.Coulomb:
            A = sqrt((sin(phi + self.varphi) * sin(phi - self.alpha)) /\
                (cos(self.varepsilon + self.varphi) * cos(self.varepsilon - self.alpha)))
            return cos(phi - self.varepsilon) ** 2 / (cos(self.varepsilon) ** 2 * \
                cos(self.varepsilon + self.varphi) * (1 + A) ** 2)

    def GetKp(self, phi, type: EarthType):
        from math import sin, cos, tan, pi, sqrt
        phi = phi * pi / 180
        if type == EarthType.Rankine:
            return tan(pi / 4 + phi / 2) ** 2
        if type == EarthType.Coulomb:
            B = sqrt((sin(phi + self.varphi) * sin(phi + self.alpha)) /\
                (cos(self.varepsilon - self.varphi) * cos(self.varepsilon - self.alpha)))
            return cos(phi + self.varepsilon) ** 2 / (cos(self.varepsilon) ** 2 * \
                cos(self.varepsilon - self.varphi) * (1 - B) ** 2)

class BoreHole(object):
    def __init__(self, soils: List[SoilMaterial], intervals:List[dict]) -> None:
        if(len(soils) != len(intervals)):
            raise Exception('The lengths of soils and intervals must be equal!')
        Soils = []
        l = len(soils)
        for i in range(l):
            dic = {
                'soil': soils[i],
                'interval': intervals[i],
            }
            Soils.append(dic)
        Soils = sorted(Soils, key=lambda i: i['interval']['top'])
        self.Soils = Soils

    def GetSoilByDeepth(self, deepth):
        for item in self.Soils:
            if item['interval']['top'] <= deepth and \
                item['interval']['bottom'] >= deepth:
                return item['soil']
        return SoilMaterial('Water', 10, 1E9, 0, 0, True)

    def GetAverageSoil(self) -> SoilMaterial:
        gamma = E = phi = c = varepsilon = varphi = alpha = intervals = thick = 0
        for item in self.Soils:
            intervals = item['interval']['bottom'] - item['interval']['top']
            gamma += item['soil'].gamma * intervals
            E += item['soil'].E * intervals
            phi += item['soil'].phi * intervals
            c += item['soil'].c * intervals
            varepsilon += item['soil'].varepsilon * intervals
            varphi += item['soil'].varphi * intervals
            alpha += item['soil'].alpha * intervals
            thick += intervals
        
        return SoilMaterial('Average soil', gamma / thick,
                            E / thick, phi / thick, c / thick,
                            True, varepsilon / thick, varphi / thick,
                            alpha / thick)


class UndergroundDiaphragmWall(object):
    """Class of the underground diaphragm wall of a foundation pit.
    """
    def __init__(self, L, h, material: ConcreteMaterial) -> None:
        """Initial the object of an underground diaphragm wall

        Args:
            L (int): Deepth of the wall, unit: m
            h (float): Thichness of the wall, unit: m
            E (int): Elastic modulus of the wall, unit: Pa
            v (float): Poisson's ratio of the wall: unit: None
            G (int): Shear modulus of the wall, unit: Pa
        """
        self.L = L
        self.h = h
        self.Material = material
        self.I = h ** 3 / 12
        self.E = material.E
        # self.EI: Stiffness of the wall
        self.EI = self.E * self.I
        # self.kar: Shear non-uniformity coefficient of Timoshenko beam
        self.kar = 5 / 6
        self.G = material.G
        self.H = 0
    
    def SetH(self, H: int):
        self.H = H

    def __str__(self) -> str:
        return """The parameters of the underground diaphragm wall are shown below:
            Length: %d m
            Thickness: %f m
            Elastic modulus: %d kPa
            Poisson's ratio: %f
        """ % (self.L, self.h, self.Material.E, self.Material.v)

class LRSide(Enum):
    LeftSide = 1
    RightSide = 2

class FoundationPit(object):
    """Class of the foundation pit, contains the parameters for calculating the instance.
    """
    def __init__(self, leftWall: UndergroundDiaphragmWall,
                 rightWall: UndergroundDiaphragmWall,
                 H1, H2, supports: List[HorizontalSupport],
                 supportCount, ds, B, D,
                 boreHole: BoreHole,
                 Palim=0.005, Pplim=0.05,
                 leftOverLoad=0, rightOverLoad=0) -> None:
        
        self.LeftWall = leftWall
        self.RightWall = rightWall
        self.H1 = H1
        self.H2 = H2
        self.ExcaveDeepth = {
            LRSide.LeftSide: H1,
            LRSide.RightSide: H2,
        }
        self.LeftWall.SetH(H1)
        self.RightWall.SetH(H2)

        self.SupportCount = supportCount
        self.Supports = supports

        self.B = B
        self.D = D

        self.BoreHole = boreHole
        self.AverageSoil = boreHole.GetAverageSoil()
        self.AverageSoil.gamma *= self.D
        
        if ds:
            self.ds = ds
        else:
            self.ds = [3 * i for i in range(supportCount)]
        
        self.UpdateLim(Palim=Palim, Pplim=Pplim)
        self.LeftOverLoad = leftOverLoad
        self.RightOverLoad = rightOverLoad

    def GetOverLoad(self, side: LRSide):
        overload = 0
        if side == LRSide.LeftSide:
            overload = self.LeftOverLoad
        if side == LRSide.RightSide:
            overload = self.RightOverLoad
        return overload * self.D

    def P0(self, z, side: LRSide):
        return self.AverageSoil.K0 * (self.AverageSoil.gamma * z + self.GetOverLoad(side))

    def Pacr(self, z, side: LRSide, type=EarthType.Rankine):
        soil = self.AverageSoil
        if type == EarthType.Rankine:
            from sympy import sqrt
            Ka = soil.Ka
            return Ka * (self.AverageSoil.gamma * z + self.GetOverLoad(side)) -\
                2 * self.AverageSoil.c * sqrt(Ka)
        if type == EarthType.Coulomb:
            from sympy import tan, atan, pi
            phi = soil.phi / 180 * pi
            H = self.ExcaveDeepth[side]
            phi = atan(tan(phi) + soil.c / soil.gamma / H)
            self.AverageSoil.phid = phi
            Ka = soil.GetKa(phi, type)
            return Ka * (self.AverageSoil.gamma * z + self.GetOverLoad(side))
        return 0

    def Ppcr(self, z, side: LRSide, type=EarthType.Rankine):
        soil = self.AverageSoil
        if type == EarthType.Rankine:
            from sympy import sqrt
            Kp = soil.Kp
            return Kp * (self.AverageSoil.gamma * z + self.GetOverLoad(side)) +\
                2 * self.AverageSoil.c * sqrt(Kp)
        if type == EarthType.Coulomb:
            from sympy import tan, atan, pi
            phi = soil.phi / 180 * pi
            H = self.ExcaveDeepth[side]
            phi = atan(tan(phi) + soil.c / soil.gamma / H)
            self.AverageSoil.phid = phi
            Kp = soil.GetKp(phi, type)
            return Kp * (self.AverageSoil.gamma * z + self.GetOverLoad(side))
        return 0

    def Pa(self, z, w, side: LRSide, alpha=0.9):
        P0 = self.P0(z, side)
        Pacr = self.Pacr(z, side)
        sa = 0
        if side ==LRSide.LeftSide:
            sa = self.PalimWl
        if side == LRSide.RightSide:
            sa = self.PalimWr
        from sympy import exp
        return P0 - (P0 - Pacr) * (w / sa) * exp(alpha * (1 - w / sa))

    def Pp(self, z, w, side: LRSide, alpha=0.9):
        H = 0
        if side == LRSide.LeftSide:
            H = self.H1
        if side == LRSide.RightSide:
            H = self.H2
        
        P0 = self.P0(z - H, side)
        Ppcr = self.Ppcr(z - H, side)
        sp = 0
        if side ==LRSide.LeftSide:
            sp = self.PplimWl
        if side == LRSide.RightSide:
            sp = self.PplimWr
        from sympy import exp
        return P0 + (Ppcr - P0) * (w / sp) * exp(alpha * (1 - w / sp))

    def UpdateLim(self, Palim=0.005, Pplim=0.05):
        self.Palim = Palim
        self.Pplim = Pplim
        self.PalimWl = self.LeftWall.L * Palim
        self.PplimWl = self.LeftWall.L * Pplim
        self.PalimWr = self.RightWall.L * Palim
        self.PplimWr = self.RightWall.L * Pplim
    
    def __str__(self) -> str:
        return """
        The foundation pit is excavated %dm on the left and %dm on the right,
        The length of left support is %dm, and the right is %dm,
        There are %d supports in it.
        """ % (self.H1, self.H2, self.LeftWall.L, self.RightWall.L, self.SupportCount)  