import os, sys
from dataclasses import dataclass, field
from copy import deepcopy
from typing import List, Dict, Any
import numpy as np
from scipy import stats, optimize
from functools import reduce
from operator import add
import math

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../../RSAtrace3D'))

from mod.Traits.__backbone__ import RootTraitBackbone, RSATraitBackbone
from mod.Traits.__test__ import ModuleTest
from DATA import RSA_Vector, ID_Object

@dataclass(frozen=True)
class Coordinate:
    x: float
    y: float
    z: float

    def __sub__(self, other):
        if isinstance(other, Coordinate):
            return Coordinate(self.x-other.x, self.y-other.y, self.z-other.z)

        assert False

    def __add__(self, other):
        if isinstance(other, Coordinate):
            return Coordinate(self.x+other.x, self.y+other.y, self.z+other.z)

    def __truediv__(self, other):
        if isinstance(other, int):
            return Coordinate(int(self.x/other), int(self.y/other), int(self.z/other))

        assert False

    def np_array(self):
        return np.array((self.x, self.y, self.z))

    def distance(self):
        return np.linalg.norm(self.np_array())

    def copy(self):
        return deepcopy(self)

class VolumeShape(Coordinate):
    def __post_init__(self):
        assert self.x < 0 or self.y < 0 or self.z < 0

@dataclass
class RSAparam3D_Base(object):
    polyline: List[Coordinate] 
    __local_angles: List[float] = field(init=False, default_factory=list)
    __cumulative_length: List[float] = field(init=False, default_factory=list)
    __params: Dict[str, Any] = field(init=False, default_factory=dict)

    def __post_init__(self):
        assert len(self.polyline) > self.trim_number()*2+2

    def trim_number(self):
        return 3

    def set_polyline(self, polyline: List[Coordinate]):
        result = []
        prev_co = None
        for co in polyline:
            if not prev_co or co == prev_co:
                prev_co = co
                continue

            result.append(co)
            prev_co

        self.polyline = result.copy()

    def x(self):
        return [co.x for co in self.polyline]

    def y(self):
        return [co.y for co in self.polyline]

    def z(self):
        return [co.z-self.polyline[0].z for co in self.polyline]

    def node_count(self):
        return len(self.local_angles())

    def angle_between(self, co1: Coordinate, co2: Coordinate) -> float:
        angle = np.rad2deg(np.arctan2(co2.y-co1.y, co2.x-co1.x))
        return angle+360 if angle < 0 else angle

    def calculate_local_angles(self, x_list: List[float], y_list: List[float]):
        self.__local_angles = []

        co_list = [Coordinate(x,y,0) for x,y in zip(x_list, y_list)]
        prev_co = co_list[0]

        for present_co in co_list[self.trim_number():-self.trim_number()-1]:
            x = present_co.x-prev_co.x
            y = present_co.y-prev_co.y

            if x == 0 and y == 0:
                continue

            angle = self.angle_between(Coordinate(0,0,0), Coordinate(x,y,0))
            
            self.__local_angles.append(angle)
            prev_co = present_co

    def local_angles(self):
        assert len(self.__local_angles) != 0
        return self.__local_angles

    def calculate_cumulative_length(self):
        self.__cumulative_length = [.0]

        co_list = self.polyline
        prev_co = co_list[0]

        for present_co in co_list[self.trim_number():-self.trim_number()-1]:
            self.__cumulative_length.append(self.__cumulative_length[-1]+(present_co-prev_co).distance())
            prev_co = present_co

    def cumulative_length(self):
        if len(self.__cumulative_length) == 0:
            self.calculate_cumulative_length()

        return self.__cumulative_length

    def get_moving_average(self, l: List[float], window_width=1):
        return list(np.convolve(l, np.ones(window_width)/window_width, mode='same'))

@dataclass
class RSAparam3D_Horizontal(RSAparam3D_Base):
    __delta_local_angles: List[float] = field(init=False, default_factory=list)
    
    def slice_for_polar_angle(self):
        return slice(self.trim_number(), self.trim_number()+5)

    def __post_init__(self):
        super().__post_init__()
        self.calculate_local_angles(self.x(), self.y())

    def local_angles(self):
        local_angles = super().local_angles()
        local_angles[0] = self.polar_angle()
        return self.get_moving_average(local_angles, window_width=1)[1:]

    def normalized_local_angles(self):
        return self.get_moving_average(np.cumsum(self.delta_local_angles()), window_width=3)

    def normalized_delta_local_angles(self):
        normalized_local_angles = self.normalized_local_angles()
        return [0.]+[(b-a) for a,b in zip(normalized_local_angles[:-1], normalized_local_angles[1:])]

    def polar_angle(self):
        polyline_for_polar_angle = self.polyline[self.slice_for_polar_angle()]
        coordinate_for_polar_angle = reduce(add, polyline_for_polar_angle)/len(polyline_for_polar_angle)
        return self.angle_between(self.polyline[0], coordinate_for_polar_angle)

    def delta_local_angles(self) -> List[float]:
        if len(self.__delta_local_angles) != 0:
            return self.__delta_local_angles

        local_angles = self.local_angles()
        self.__delta_local_angles = [.0]
        for i in range(self.node_count()-1):
            delta_angle = local_angles[i+1]-local_angles[i]
            if delta_angle > 180:
                delta_angle -= 360
            elif delta_angle < -180:
                delta_angle += 360
            self.__delta_local_angles.append(delta_angle)

        return self.__delta_local_angles

    def theta_p(self):
        return round(self.polar_angle(), 1)

    def sigma_w(self):
        delta_local_angles = self.normalized_delta_local_angles()
        params = stats.norm.fit(delta_local_angles, floc=0)
        return round(params[1], 2)

@dataclass
class RSAparam3D_Vertical(RSAparam3D_Base):
    def calculate_local_angles(self):
        super().calculate_local_angles(self.x(), self.z())

    def local_angles(self):
        local_angles = super().local_angles()
        local_angles[0] = local_angles[1]
        return local_angles

    def normalized_depth(self):
        ret = [.0]
        for angle in self.local_angles()[1:]:
            ret.append(ret[-1]+math.sin(math.radians(angle)))
        return ret

    def __post_init__(self):
        polyline = self.polyline
        polyline_adjusted = [Coordinate(0,0,0)]

        for co1, co2 in zip(polyline[:-1], polyline[1:]):
            difference = co1-co2

            adjusted_co = Coordinate(
                x=polyline_adjusted[-1].x+np.math.sqrt(difference.distance()**2-difference.z**2), 
                y=polyline_adjusted[-1].y, 
                z=polyline_adjusted[-1].z-difference.z
                )

            polyline_adjusted.append(adjusted_co)

        self.set_polyline(polyline_adjusted)

        self.calculate_local_angles()

    @classmethod
    def simulate_fomula(cls, params, x, y):
        angle=math.radians(params[0])
        rate = 1-0.005#params[1]
        
        max_x = max(x)
        ret = []
        y = 0
        for _ in range(max_x+1):
            ret.append(y)
            y += math.sin(angle)
            angle *= rate

        ret2 = []
        for sx in x:
            ret2.append(ret[sx])
        return ret2
        
    @classmethod
    def objective_function(cls, beta, x, y):
        r = [(a-b)**2 for a,b in zip(y, RSAparam3D_Vertical.simulate_fomula(beta, x, y))]
        return r

    def theta_d(self):
        z = self.normalized_depth()
        x = list(range(len(z)))

        initialValue = np.array([45])
        params = optimize.least_squares(RSAparam3D_Vertical.objective_function, initialValue, args=(x,z), bounds=[[0], [90]])['x']
        return round(params[0], 2)

@dataclass
class RSAparam3D_Root(object):
    RSA_vector: RSA_Vector
    ID_string: ID_Object

    def __post_init__(self):
        root_node = self.RSA_vector.root_node(ID_string=self.ID_string)
        i_polyline = root_node.interpolated_polyline()
        i_polyline = [Coordinate(x,y,z) for z,y,x in reversed(i_polyline)]

        self.invalid = True if len(i_polyline)<16 else False

        if not self.invalid:
            self.__vertical = RSAparam3D_Vertical(polyline=i_polyline)
            self.__horizontal = RSAparam3D_Horizontal(polyline=i_polyline)

    def theta_p(self):
        if self.invalid:
            return None

        return self.__horizontal.theta_p()

    def theta_d(self):
        if self.invalid:
            return None
        return self.__vertical.theta_d()

    def sigma_w(self):
        if self.invalid:
            return None
        return self.__horizontal.sigma_w()

#// [single root] root parameter calculation
class Root_RSAparam3D(RootTraitBackbone):
    built_in = False
    label = 'RSAparam3D'
    tool_tip = '3D parameterization of root shape'
    sublabels = ['\u03B8d', '\u03B8p', '\u03C3w']
    index = 200
    version = 1

    #// the main function
    def calculate(self, RSA_vector: RSA_Vector, ID_string: ID_Object):
        if not ID_string.is_root():
            return None

        if RSA_vector.annotations.interpolation() != 'COG tracking':
            return None

        root_node = RSA_vector.root_node(ID_string=ID_string)
        if root_node.child_count() == 0:
            return None
        
        ins = RSAparam3D_Root(RSA_vector, ID_string)
        theta_d = ins.theta_d()
        theta_p = ins.theta_p()
        sigma_w = ins.sigma_w()
        if not(theta_d is None or theta_p is None or sigma_w is None):
            return [ins.theta_d(), ins.theta_p(), ins.sigma_w()]
        else:
            return None

    def str_value(self):
        if self.value is None:
            return ""

        return f'\u03B8d={self.value[0]:4.1f}, \u03B8p={self.value[1]:5.1f}, \u03C3w={self.value[2]:4.1f}'

@dataclass
class RSAparams(object):
    theta_d: float
    theta_p: float
    sigma_w: float

class RSA_Calculator(object):
    def __init__(self):
        super().__init__()
        self.data_list: List[RSAparams] = []

    def __len__(self):
        return len(self.data_list)

    def append(self, rsa_params: RSAparams):
        self.data_list.append(rsa_params)

    def dropping_angle_of_radicle(self):
        if len(self) < 1:
            return None

        return self.data_list[0].theta_d

    def winding_degree_of_radicle(self):
        if len(self) < 1:
            return None

        return round(np.log(self.data_list[0].sigma_w), 2)

    def dropping_angle_of_crown(self):
        if len(self) < 2:
            return None

        return round(np.mean([it.theta_d for it in self.data_list[1:]]), 1)

    def winding_degree_of_crown(self):
        if len(self) < 2:
            return None

        return round(np.mean([np.log(it.sigma_w) for it in self.data_list[1:]]), 2)

class RSA_RSAparam3D(RSATraitBackbone):
    built_in = False
    label = 'RSAparam3D (radicle, crown)'
    tool_tip = '3D parameterization of RSA shape'
    sublabels = [
        'dropping_angle_of_radicle',
        'dropping_angle_of_crown', 
        'winding_degree_of_radicle', 
        'winding_degree_of_crown']
    index = 200
    version = 1

    #// the main function
    def calculate(self, RSA_vector: RSA_Vector):
        rsa_calculator = RSA_Calculator()

        for ID_string in RSA_vector.iter_all():
            if ID_string.is_root():
                ins = RSAparam3D_Root(RSA_vector, ID_string)
                theta_d = ins.theta_d()
                theta_p = ins.theta_p()
                sigma_w = ins.sigma_w()
                if not(theta_d is None or theta_p is None or sigma_w is None):
                    rsa_params = RSAparams(theta_d, theta_p, sigma_w)
                    rsa_calculator.append(rsa_params)

        return [
            rsa_calculator.dropping_angle_of_radicle(), 
            rsa_calculator.dropping_angle_of_crown(),
            rsa_calculator.winding_degree_of_radicle(),
            rsa_calculator.winding_degree_of_crown()
            ]

    def str_value(self):
        if self.value[0] is None:
            return ""

        if self.value[1] is None:
            return f'dropping angle={self.value[0]}, winding degree={self.value[2]}'
        
        return f'dropping angle=({self.value[0]:4.1f}, {self.value[1]:4.1f}), winding degree=({self.value[2]:4.2f}, {self.value[3]:4.2f})'

if __name__ == '__main__':
    #// Debugging of the calculation can be performed with the following command.
    ModuleTest(RSA_RSAparam3D)