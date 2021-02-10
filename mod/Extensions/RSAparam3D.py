from copy import deepcopy
from dataclasses import dataclass, field
from itertools import accumulate
import logging
import numpy as np
import operator
import os
import matplotlib.pyplot as plt
from scipy import stats
import sys
from typing import List, Tuple, Union

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../../RSAtrace3D'))
    logging.basicConfig(level=logging.INFO)

from PyQt5.QtWidgets import QCheckBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QDoubleValidator
from mod.Extensions.__backbone__ import ExtensionBackbone
from GUI import QtMain
from DATA import RSA_Vector

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

@dataclass
class CoordicateVector(object):
    coordinates: List[Coordinate] = field(default_factory=list)

    def append(self, coordicate: Coordinate):
        self.coordinates.append(coordicate)

    def x(self):
        return [co.x for co in self.coordinates]

    def y(self):
        return [co.y for co in self.coordinates]

    def z(self):
        return [co.z for co in self.coordinates]

    def xyz(self):
        return [self.x(), self.y(), self.z()]

    def xyz_mm_unit(self, resolution: float):
        return (np.asarray(self.xyz())*resolution).tolist()

class VolumeShape(Coordinate):
    def __post_init__(self):
        assert self.x < 0 or self.y < 0 or self.z < 0

@dataclass(frozen=True)
class RSAparam3D_Params(object):
    root_number: int
    dropping_angle_crown: float
    dropping_angle_sd_crown: float
    winding_degree_crown: float
    winding_degree_sd_crown: float
    dropping_angle_radicle: float
    winding_degree_radicle: float
    seed: int
    resolution: float = 0.3
    iteration: int = 10000
    depth: int = 860
    radius: int = 300

@dataclass
class RSAparam3D_Root(object):
    polar_angle: float
    dropping_angle: float
    log_mu_w: float
    depth: int
    radius: int
    max_iteration: int
    size: int = 4

    def __post_init__(self):
        dropping_angle = self.dropping_angle

        frozen_norm_winding = stats.norm.freeze(loc=0, scale=np.exp(self.log_mu_w))
        polar_angles  = [self.polar_angle+a for a in accumulate(frozen_norm_winding.rvs(self.max_iteration))]
        dropping_angles = [self.dropping_angle*factor for factor in accumulate([1-0.005]*self.max_iteration, operator.mul)]

        self.vector = CoordicateVector()
        z, y, x = (.0,.0,.0)
        for i in range(self.max_iteration):
            self.vector.append(Coordinate(x,y,z))

            if i == 0:
                self.vector.append(Coordinate(x,y,z))
            else:
                x += np.cos(np.deg2rad(polar_angles[i]))*self.size*np.cos(np.deg2rad(dropping_angles[i]))
                y += np.sin(np.deg2rad(polar_angles[i]))*self.size*np.cos(np.deg2rad(dropping_angles[i]))
                
                z += np.sin(np.deg2rad(dropping_angles[i]))*self.size

                if z > self.depth or np.sqrt(x**2+y**2) > self.radius:
                    break

        #print(self.vector)

        #self.vector = [(int(z), int(y), int(x)) for z,y,x in self.vector]

        #x = [x for z, y, x in self.vector]
        #y = [y for z, y, x in self.vector]



@dataclass
class RSAparam3D_RSA(object):
    parameters: RSAparam3D_Params

    def __post_init__(self):
        self.roots: List[RSAparam3D_Root] = []

    def get_shape(self):
        return VolumeShape(self.parameters.radius*2-1, self.parameters.radius*2-1, self.parameters.depth)

    def is_calculated(self):
        return len(self.roots) != 0

    def calculate(self):
        np.random.seed(self.parameters.seed)

        #// calculate parameters of each root
        self.list_of_polar_angles = list(np.random.uniform(0,360, self.parameters.root_number))
        self.list_of_log_mu_w = list(np.random.normal(self.parameters.winding_degree_crown, self.parameters.winding_degree_sd_crown, self.parameters.root_number).clip(0))
        self.list_of_log_mu_w[0] = self.parameters.winding_degree_radicle
        self.list_of_dropping_angles = list(np.random.normal(self.parameters.dropping_angle_crown, self.parameters.dropping_angle_sd_crown, self.parameters.root_number).clip(0,90))
        self.list_of_dropping_angles[0] = self.parameters.dropping_angle_radicle

        for i in range(self.parameters.root_number):
            self.roots.append(RSAparam3D_Root(
                polar_angle=self.list_of_polar_angles[i], 
                dropping_angle=self.list_of_dropping_angles[i],
                log_mu_w=self.list_of_log_mu_w[i],
                depth=self.parameters.depth,
                radius=self.parameters.radius, 
                max_iteration=self.parameters.iteration
            ))

    def simulate(self):
        if not self.is_calculated():
            self.calculate()
        
        figure = plt.figure()
        ax = plt.axes(projection="3d")
        ax.set_xlim([-self.parameters.radius*self.parameters.resolution, self.parameters.radius*self.parameters.resolution])
        ax.set_ylim([-self.parameters.radius*self.parameters.resolution, self.parameters.radius*self.parameters.resolution])
        ax.set_zlim([self.parameters.depth*self.parameters.resolution, 0])
        ax.set_box_aspect((1,1,1))

        for root in self.roots:
            x, y, z = root.vector.xyz_mm_unit(resolution=self.parameters.resolution)
            ax.plot3D(x, y, z)

            

        plt.show()

    def make_rinfo(self, volume_name:str):
        if not self.is_calculated():
            self.calculate()

        base_loc = Coordinate(self.parameters.radius//2, self.parameters.radius//2, 0)

        rinfo = RSA_Vector()
        rinfo.annotations.set_resolution(0.3)
        rinfo.annotations.set_volume_shape((base_loc.z, base_loc.y, base_loc.x))
        rinfo.annotations.set_volume_name(volume_name)
        rinfo.annotations.set_interpolation('RSAparam3D simulator')

        annotations={'coordinate': [base_loc.z, base_loc.y, base_loc.x]}
        base_ID_string = rinfo.append_base(annotations=annotations)

        for i, root in enumerate(self.roots):
            polyline = [[int(z),int(x+base_loc.x),int(y+base_loc.y)] for x,y,z in zip(root.vector.x(), root.vector.y(), root.vector.z())][::-1]
            annotations={'polyline': polyline}
            base_node = rinfo.base_node(ID_string=base_ID_string)
            root_ID_string = base_node.append(annotations=annotations)
            root_node = rinfo.root_node(ID_string=root_ID_string)
            annotations={'coordinate': polyline[0]}

            root_node.append(annotations=annotations, interpolation=False)

        return rinfo

@dataclass
class RSAparam3D_Simulator(object):
    parameters: RSAparam3D_Params

    def simulate(self):
        rsa = RSAparam3D_RSA(self.parameters)
        rsa.simulate()

    def make_rinfo(self, volume_name:str):
        rsa = RSAparam3D_RSA(self.parameters)
        return rsa.make_rinfo(volume_name=volume_name)

class LineEditLimited(QLineEdit):
    def __init__(self, *args, limits: Tuple[Union[int, float], Union[int, float]]=None, **kwargs):
        super().__init__(*args, **kwargs)
        validator = QDoubleValidator()
        validator.setBottom(0.)
        self.setValidator(validator)

        self.limits = limits or (0, 90)

        self.textChanged.connect(self.on_text_changed)

    def postFocusOutEvent(self, ev):
        return

    def focusOutEvent(self, ev):
        super().focusOutEvent(ev)
        if self.text().endswith('.'):
            self.setText(f'{self.text()[:-1]}')

        if self.text() == '':
            self.setText(f'{self.limits[0]}')

        self.postFocusOutEvent(ev)
        
    def on_text_changed(self, text:str):
        if 'e' in text:
            text = text.replace('e', '')
            self.setText(text)

        if text == '' or text == '.':
            return
            
        if type(self.limits[1]) is float and float(text) > self.limits[1]:
                self.setText(f'{self.limits[1]}')
        if type(self.limits[1]) is int and np.ceil(float(text)) > self.limits[1]:
                self.setText(f'{self.limits[1]}')

        if type(self.limits[0]) is float and float(text) < self.limits[0]:
                self.setText(f'{self.limits[0]}')
        if type(self.limits[0]) is int and np.floor(float(text)) < self.limits[0]:
                self.setText(f'{self.limits[0]}')

    def keyPressEvent(self, ev):
        if ev.key() == Qt.Key_Return and not ev.isAutoRepeat():
            ev.accept()
            self.clearFocus()

        return super().keyPressEvent(ev)

class LineEdit_RootNumber(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(1,100), *args, **kwargs)

    def on_text_changed(self, text:str):
        if '.' in text:
            text = text.replace('.', '')
            self.setText(text)

        super().on_text_changed(text)

class LineEdit_Seed(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(0,1000), *args, **kwargs)

    def on_text_changed(self, text:str):
        if '.' in text:
            text = text.replace('.', '')
            self.setText(text)

        super().on_text_changed(text)

class LineEdit_DroppingMean(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(0,80), *args, **kwargs)

class LineEdit_DroppingSD(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(0,30), *args, **kwargs)

class LineEdit_WindingMean(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(0.1,3), *args, **kwargs)

class LineEdit_WindingSD(LineEditLimited):
    def __init__(self, *args, **kwargs):
        super().__init__(limits=(0.1,1), *args, **kwargs)

class Group_Dropping(object):
    def __init__(self):
        self.items = (
            QLabel('<b>Crown root:</b>'),
            QLabel('dropping angle'), 
            LineEdit_DroppingMean('50'), 
            QLabel('\u00b1'),
            LineEdit_DroppingSD('40'), 
            QCheckBox('auto')
        )

        self.items[5].stateChanged.connect(self.on_checkbox_changed)
        self.items[5].setChecked(True)

    def on_checkbox_changed(self, state):
        if state == 2:
            self.items[4].setReadOnly(True)

            item: QLineEdit = self.items[2]
            v =  float(item.text())*0.144+6.518
            
            palette = QPalette()
            palette.setColor(QPalette.Base, Qt.gray)
            palette.setColor(QPalette.Text, Qt.darkGray)
            self.items[4].setPalette(palette)

            self.items[4].setText(f'{v:.2f}')
        else:
            self.items[4].setReadOnly(False)
            
            palette = QPalette()
            palette.setColor(QPalette.Base, Qt.white)
            palette.setColor(QPalette.Text, Qt.black)
            self.items[4].setPalette(palette)

class Group_RootNumber(object):
    def __init__(self):
        self.items = (
            QLabel('<b>General:</b>'),
            QLabel('root number'), 
            LineEdit_RootNumber('25'), 
            QLabel(''),
            QLabel(''),
            QLabel('')
        )

class Group_Seed(object):
    def __init__(self):
        self.items = (
            QLabel(''),
            QLabel('random seed'), 
            LineEdit_Seed('1'), 
            QLabel(''),
            QLabel(''),
            QLabel('')
        )

class Group_Winding(object):
    def __init__(self):
        self.items = (
            QLabel(''),
            QLabel('winding degree'), 
            LineEdit_WindingMean('1.89'), 
            QLabel('\u00b1'),
            LineEdit_WindingSD('3.1'), 
            QCheckBox('auto')
        )

        self.items[5].stateChanged.connect(self.on_checkbox_changed)
        self.items[5].setChecked(True)

    def on_checkbox_changed(self, state):
        if state == 2:
            self.items[4].setReadOnly(True)

            item: QLineEdit = self.items[2]
            v =  float(item.text())*0.348-0.302
            
            palette = QPalette()
            palette.setColor(QPalette.Base, Qt.gray)
            palette.setColor(QPalette.Text, Qt.darkGray)
            self.items[4].setPalette(palette)

            self.items[4].setText(f'{v:.2f}')
        else:
            self.items[4].setReadOnly(False)
            
            palette = QPalette()
            palette.setColor(QPalette.Base, Qt.white)
            palette.setColor(QPalette.Text, Qt.black)
            self.items[4].setPalette(palette)

class Group_Radicle_Dropping(object):
    def __init__(self):
        self.items = (
            QLabel('<b>Radicle:</b>'),
            QLabel('dropping angle'), 
            DroppingAngleLineEdit('50'), 
            QLabel(''),
            QLabel(''), 
            QLabel('')
        )

class Group_Radicle_Winding(object):
    def __init__(self):
        self.items = (
            QCheckBox('auto'),
            QLabel('winding degree'), 
            LineEdit_WindingMean('3.1'), 
            QLabel(''),
            QLabel(''), 
            QLabel('')
        )

class ParameterPanel(QWidget):
    def simulate(self):
        parameters = RSAparam3D_Params(**self.get_params())
        simulator = RSAparam3D_Simulator(parameters)
        simulator.simulate()

    def get_params(self):
        return {
            'root_number': int(self.g_root_number.items[2].text()), 
            'seed': int(self.g_seed.items[2].text()), 
            'dropping_angle_crown': float(self.g_crown_dropping.items[2].text()),
            'dropping_angle_sd_crown': float(self.g_crown_dropping.items[4].text()),
            'winding_degree_crown': float(self.g_crown_winding.items[2].text()),
            'winding_degree_sd_crown': float(self.g_crown_winding.items[4].text()),
            'dropping_angle_radicle': float(self.g_radicle_dropping.items[2].text()),
            'winding_degree_radicle': float(self.g_radicle_winding.items[2].text()),
        }

    def on_rinfo_button_clicked(self):
        save_file_name, _ = QFileDialog.getSaveFileName(self, None, os.path.expanduser('~'), 'rinfo file (*.rinfo)')
        if save_file_name == '':
            return
            
        if not save_file_name.lower().endswith('.rinfo'):
            save_file_name += '.rinfo'

        parameters = RSAparam3D_Params(**self.get_params())
        simulator = RSAparam3D_Simulator(parameters)
        rinfo = simulator.make_rinfo(os.path.basename(save_file_name))
        rinfo.save(rinfo_file_name=save_file_name)

    def __init__(self, parent):
        super().__init__(**{'parent': parent})

        self.horizontal_box = QHBoxLayout()

        self.sim_button = QPushButton("Simulate")
        self.sim_button.clicked.connect(self.simulate)
        self.rinfo_button = QPushButton("Save as .rinfo")
        self.rinfo_button.clicked.connect(self.on_rinfo_button_clicked)

        self.horizontal_box.addStretch(1)
        self.horizontal_box.addWidget(self.sim_button)
        self.horizontal_box.addWidget(self.rinfo_button)

        self.vertical_box = QVBoxLayout()

        self.g_root_number = Group_RootNumber()
        self.g_seed = Group_Seed()
        self.g_crown_dropping = Group_Dropping()
        self.g_crown_winding = Group_Winding()
        self.g_radicle_dropping =  Group_Radicle_Dropping()
        self.g_radicle_winding = Group_Radicle_Winding()

        def on_checkbox_changed(state):
            if state == 2:
                self.g_radicle_dropping.items[2].setReadOnly(True)
                self.g_radicle_winding.items[2].setReadOnly(True)

                item1 = self.g_crown_dropping.items[2]
                item2 = self.g_crown_winding.items[2]
                v1 =  float(item1.text())*0.241+49.6
                v2 = float(item2.text())*0.969+1.13
                
                palette = QPalette()
                palette.setColor(QPalette.Base, Qt.gray)
                palette.setColor(QPalette.Text, Qt.darkGray)
                self.g_radicle_dropping.items[2].setPalette(palette)
                self.g_radicle_winding.items[2].setPalette(palette)

                self.g_radicle_dropping.items[2].setText(f'{v1:.2f}')
                self.g_radicle_winding.items[2].setText(f'{v2:.2f}')
            else:
                self.g_radicle_dropping.items[2].setReadOnly(False)
                self.g_radicle_winding.items[2].setReadOnly(False)
                
                palette = QPalette()
                palette.setColor(QPalette.Base, Qt.white)
                palette.setColor(QPalette.Text, Qt.black)
                self.g_radicle_dropping.items[2].setPalette(palette)
                self.g_radicle_winding.items[2].setPalette(palette)

        self.g_radicle_winding.items[0].stateChanged.connect(on_checkbox_changed)
        self.g_radicle_winding.items[0].setChecked(True)

        def postFocusOutEvent(ev):
            if self.g_crown_dropping.items[5].isChecked():
                self.g_crown_dropping.on_checkbox_changed(state=2)
            if self.g_crown_winding.items[5].isChecked():
                self.g_crown_winding.on_checkbox_changed(state=2)
            if self.g_radicle_winding.items[0].isChecked():
                on_checkbox_changed(state=2)

        setattr(self.g_crown_dropping.items[2], 'postFocusOutEvent', postFocusOutEvent)
        setattr(self.g_crown_winding.items[2], 'postFocusOutEvent', postFocusOutEvent)
    

        g_list = [
            self.g_root_number, 
            self.g_seed,
            self.g_crown_dropping, 
            self.g_crown_winding,
            self.g_radicle_dropping,
            self.g_radicle_winding
        ]

        vertical_box_list = [QVBoxLayout() for _ in range(6)]

        for g in g_list:
            for vertical_box, it in zip(vertical_box_list, g.items):
                vertical_box.addWidget(it)

        hbox = QHBoxLayout()
        for vertical_box in vertical_box_list:
            hbox.addLayout(vertical_box)

        hbox.addStretch(1)

        self.vertical_box.addLayout(hbox)

        self.vertical_box.addLayout(self.horizontal_box)
        self.vertical_box.addStretch(1)

        self.horizontal_box_final = QHBoxLayout()
        self.horizontal_box_final.addLayout(self.vertical_box)
        self.horizontal_box_final.addStretch(1)

        self.setLayout(self.horizontal_box_final)

class ExtensionRSAparam3D(QMainWindow, ExtensionBackbone):
    label = 'RSAparam3D simulator'
    status_tip = 'compute RSA with parameters'
    index = 200
    version = 1

    def __init__(self, parent: QtMain = None):
        super().__init__(parent=parent)

        self.setWindowTitle('RSAparam3D simulator')
        self.__parent = parent

        self.param_panel = ParameterPanel(self)

        self.setCentralWidget(self.param_panel)

class DroppingAngleLineEdit(QLineEdit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setToolTip('0 to 90')
        validator = QDoubleValidator()
        validator.setBottom(0.)
        self.setValidator(validator)

        self.limits = (0, 90)

        self.textChanged.connect(self.on_text_changed)

    def focusOutEvent(self, ev):
        super().focusOutEvent(ev)
        if self.text().endswith('.'):
            self.setText(f'{self.text()[:-1]}')

        if self.text() == '':
            self.setText(f'{self.limits[0]}')
        
    def on_text_changed(self, text:str):
        if 'e' in text:
            text = text.replace('e', '')
            self.setText(text)

        if text == '' or text == '.':
            return
            
        if type(self.limits[1]) is float and float(text) > self.limits[1]:
                self.setText(f'{self.limits[1]}')
        if type(self.limits[1]) is int and np.ceil(float(text)) > 90:
                self.setText(f'{self.limits[1]}')

    def keyPressEvent(self, ev):
        if ev.key() == Qt.Key_Return and not ev.isAutoRepeat():
            ev.accept()
            self.clearFocus()

        return super().keyPressEvent(ev)

if __name__ == '__main__':
    from PyQt5.QtWidgets import QApplication
    app = QApplication([])
    ins = ExtensionRSAparam3D()
    ins.show()
    app.exec_()