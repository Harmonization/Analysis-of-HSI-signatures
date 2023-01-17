from PyQt5 import QtWidgets as QW, QtGui, QtCore
from Form import Ui_MainWindow  # импорт нашего сгенерированного файла

import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.colors as mcolors, cv2
import re, os, sys, pysptools.util as util, pylab, matplotlib

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NTool
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
matplotlib.use('Qt5Agg')
wave = np.load('wave.npy')[:204]
COLOR = 'nipy_spectral'

class mywindow(QW.QMainWindow):
    def __init__(self):
        super(mywindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.ui.button_open.clicked.connect(self.Open_HSI)
        self.showMaximized()
        self.Open_HSI()
    
    def Open_HSI(self):
        path = QW.QFileDialog.getOpenFileName(self, "Open File", '', '*.hdr')[0]
        #path = 'C:\\Users\\inter\\0. Изучение\\Experiment_1\\Data\\1_45\\3\\1.hdr'
        self.hsi = np.rot90(util.load_ENVI_file(path)[0], k = 3)
        #self.hsi[self.hsi > 1] = 1
        self.rgb = self.hsi[:, :, (70, 53, 19)]
        self.signatures = []
        
        # Visualisation
        fig_hsi, ax_hsi = plt.subplots(); fig, ax = plt.subplots()
        fig_hsi.subplots_adjust(left=0, right=1, bottom=0, top=1)
        fig.subplots_adjust(left=.05, right=.95, bottom=.05, top=.95)
        ax_hsi.imshow(self.rgb); ax_hsi.axis('off')
        
        canvas_hsi = FigureCanvasQTAgg(fig_hsi); canvas = FigureCanvasQTAgg(fig)
        toolbar_hsi = NTool(canvas_hsi, self); toolbar = NTool(canvas, self)
        
        # Buttons
        next_button = QW.QPushButton('Выбрать ROI')
        next_button.clicked.connect(self.Edit_ROI)
        
        def Clear():
            ax_hsi.clear(); ax.clear()#; ax.set_ylim([-.05, 1.05])
            ax_hsi.imshow(self.rgb); ax_hsi.axis('off')
            canvas_hsi.draw(); canvas.draw()
            self.signatures = []
        
        clear_button = QW.QPushButton('Сбросить')
        clear_button.clicked.connect(Clear)
        
        def Save():
            # Сохранить спектры в таблицу
            if self.signatures:
                df = pd.DataFrame(self.signatures, columns=wave)
                save_path = ''
                save_path = ''.join(QW.QFileDialog.getSaveFileName(self, "Save File", '', '.xlsx'))
                if save_path:
                    df.to_excel(save_path)
            
        
        save_button = QW.QPushButton('Сохранить')
        save_button.clicked.connect(Save)
        
        # Комбинироанный виджет
        hbox = QW.QHBoxLayout()
        vbox_1, vbox_2 = QW.QVBoxLayout(), QW.QVBoxLayout()
        vbox_1.addWidget(canvas_hsi); vbox_1.addWidget(toolbar_hsi)
        vbox_2.addWidget(next_button); vbox_2.addWidget(clear_button)
        vbox_2.addWidget(save_button); vbox_2.addWidget(self.ui.button_open)
        vbox_2.addWidget(canvas); vbox_2.addWidget(toolbar)
        hbox.addLayout(vbox_1); hbox.addLayout(vbox_2)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QW.QWidget()
        widget.setLayout(hbox)
        self.setCentralWidget(widget)
        self.show()
        
        # Interact
        def onclick(event):
            if str(event.button) == 'MouseButton.RIGHT':
                x, y = self.coord(event) # координаты пикселя
                ax_hsi.scatter(x=x, y=y) # Точка на изображении
                ax.plot(wave, self.hsi[y, x, :], linewidth=2) # Спектр на графике
                canvas_hsi.draw(); canvas.draw()
                
                self.signatures.append(self.hsi[y, x, :])
            
        fig_hsi.canvas.mpl_connect('button_press_event', onclick)
    
    def coord(self, txt):
        t = re.findall("\d+\.\d+", str(txt))
        return round(float(t[0])), round(float(t[1]))
        
    def Edit_ROI(self):
        self.mpl_rect = []; self.matr_a = 3; self.matr_b = 99; self.band_1 = 18; self.band_2 = 54; self.Xy = {'mpl_rects': []}
        self.cbar = 0
        # Отделение растения от фона
        self.w = self.hsi[:, :, self.band_2] - self.hsi[:, :, self.band_1]
        self.w /= (self.w.max() - self.w.min())
        
        fig_l, ax_l = plt.subplots(); fig_r, ax_mx = plt.subplots(); ax_ms = fig_r.add_axes([.18, .7, .4, .2])
        fig_l.subplots_adjust(left=0, right=1, bottom=0, top=1)
        fig_r.subplots_adjust(left=.05, right=.95, bottom=.05, top=.95)
        
        canvas_l = FigureCanvasQTAgg(fig_l); canvas_r = FigureCanvasQTAgg(fig_r)
        toolbar_l = NTool(canvas_l, self); toolbar_r = NTool(canvas_r, self)
        
        # B U T T O N S
        def Save():
            # Сохранить спектры в таблицу
            save_path = ''
            save_path = QW.QFileDialog.getSaveFileName(self, "Save File", '', '.xlsx')
            if save_path:
                ms_lst = []
                for xy_num, xy_lst in self.Xy.items():
                    if type(xy_num) == int and len(xy_lst) == 4:
                        x1, y1, x2, y2 = xy_lst
                        mask = self.mask['rect'][y1: y2, x1: x2]
                        ms = MSpectre(x1, y1, x2, y2, mask)
                        ms_lst.append(list(ms))
                        mx = pd.DataFrame(MMatrix(ms, vis=False), index=wave[self.matr_a: self.matr_b + 1], 
                                          columns=wave[self.matr_a: self.matr_b + 1])
                        mx.to_excel(save_path[0] + f'_matrix_{xy_num}' + save_path[1])
                
                ms = pd.DataFrame(ms_lst, columns=wave)
                ms.to_excel(save_path[0] + '_spectre' + save_path[1])
        
        save_button = QW.QPushButton('Сохранить')
        save_button.clicked.connect(Save)
        
        
        rb_1 = QW.QRadioButton('Изображение целиком'); rb_2 = QW.QRadioButton('Только ROI'); rb_1.setChecked(True)
        
        def Text_Changed():
            a, b = int(text_1.text()) - 1, int(text_2.text()) - 1
            if a < 0: a = 0
            if a > 202: a = 202
            if b < 1: b = 1
            if b > 203: b = 203
            if a >= b: 
                a, b = self.matr_a, self.matr_b
            text_1.setText(str(a + 1)); text_2.setText(str(b + 1))
            self.matr_a, self.matr_b = a, b
            self.mmatrix = np.zeros((self.matr_a, self.matr_b + 1))
            Draw_Right()
            
        text_1 = QW.QLineEdit(str(self.matr_a + 1)); text_1.returnPressed.connect(Text_Changed); text_1.setFixedSize(35, 35)
        text_2 = QW.QLineEdit(str(self.matr_b + 1)); text_2.returnPressed.connect(Text_Changed); text_2.setFixedSize(35, 35)
        
        def Text_W_Changed():
            a, b = int(text_w_1.text()) - 1, int(text_w_2.text()) - 1
            if a < 0: a = 0
            if a > 203: a = 203
            if b < 0: b = 0
            if b > 203: b = 203
            text_w_1.setText(str(a + 1)); text_w_2.setText(str(b + 1))
            self.band_1, self.band_2 = a, b
            
            self.w = self.hsi[:, :, self.band_2] - self.hsi[:, :, self.band_1]
            self.w /= (self.w.max() - self.w.min())
            #thr_slider.setMinimum = int(self.w.min() * 100 - 1); thr_slider.setMaximum = int(self.w.max() * 100 + 1)
            #thr_slider.setValue(int(np.percentile(self.w, 85) * 100))
            Filter()
        
        text_w_1 = QW.QLineEdit(str(self.band_1 + 1)); text_w_1.returnPressed.connect(Text_W_Changed); text_w_1.setFixedSize(35, 35)
        text_w_2 = QW.QLineEdit(str(self.band_2 + 1)); text_w_2.returnPressed.connect(Text_W_Changed); text_w_2.setFixedSize(35, 35)
        
        def MMatrix(sign, vis=True):
            # Calculate
            s = sign[self.matr_a: self.matr_b + 1]
            mx = np.tril(s[:, np.newaxis] - s).T
            if not vis: return mx
            
            # Vis
            ax_mx.clear(); ax_mx.imshow(mx, cmap=COLOR, aspect='equal', vmin=-1, vmax=1, origin="lower")
            ax_indx = np.linspace(self.matr_a, self.matr_b, 6).astype(int)
            ticks = ax_indx - self.matr_a; ax_mx.set_xticks(ticks); ax_mx.set_yticks(ticks)
            ax_mx.set_xticklabels(wave[ax_indx]); ax_mx.set_yticklabels(wave[ax_indx])
            if not self.cbar: self.cbar = fig_r.colorbar(ScalarMappable(norm=mcolors.Normalize(vmin=-1, vmax=1), 
                                                                    cmap=COLOR), ax=ax_mx)
        def MSpectre(x1, y1, x2, y2, mask):
            # Calculate
            ms = self.hsi[y1: y2, x1: x2][mask].mean(axis=0)
            #Vis
            ax_ms.clear(); ax_ms.plot(wave, ms, color='black', linewidth=2)#; ax_ms.set_ylim([-.05, 1.05])
            ax_ms.set_title('Средний спектр', fontsize=11)
            return ms
            
        def Draw_Right():
            # Ms + Mx
            x1, y1, x2, y2 = self.Xy[combo_box.currentData()]
            mask = self.mask['rect'][y1: y2, x1: x2]
            if np.any(mask):
                MMatrix(MSpectre(x1, y1, x2, y2, mask))
            else:
                ax_ms.clear(); ax_mx.clear(); ax_ms.axis('off'); ax_mx.axis('off')
            canvas_r.draw()
            
        combo_box = QW.QComboBox(); combo_box.activated.connect(Draw_Right)    
        
        def Calculate_mask(T, kick=1):
            # Calculate mask, kick=0 если учитывать маску всего изображения
            mask_full = np.zeros(self.rgb.shape[:2], dtype=bool)
            mask_full[self.w >= T] = True
            
            mask_rect = np.zeros(self.rgb.shape[:2], dtype=bool)
            for xy_num, xy_lst in self.Xy.items():
                if type(xy_num) == int and len(xy_lst) == 4:
                    if kick == 1 and xy_num == 0: continue
                    x1, y1, x2, y2 = xy_lst
                    mask_rect[y1: y2, x1: x2][mask_full[y1: y2, x1: x2]] = True
                        
            self.mask = {'full': mask_full, 'rect': mask_rect}
            
        def Draw_Left():
            key = 'full' if rb_1.isChecked() else 'rect'
            img = self.rgb.copy(); img[self.mask[key] == 0] = 0
            if self.Xy['ax_l']: self.Xy['ax_l'].remove()
            self.Xy['ax_l'] = ax_l.imshow(img); ax_l.axis('off'); canvas_l.draw()
        
        def Filter():
            # Событие при изменении ползунка
            T = thr_slider.value() * .01
            Calculate_mask(T, kick=0 if combo_box.currentData() == 0 else 1) # Вычисление маски
            Draw_Left() # Отрисовка rgb
            Draw_Right() # Отрисовка правой части
            label_slider.setText(f'{np.around(T, 2)}')
        
        def Change_Thr():
            label_slider.setText(f'{np.around(thr_slider.value() * .01, 2)}')
        
        
        thr_slider = QW.QSlider(minimum=int(self.w.min() * 100 - 1), maximum=int(self.w.max() * 100 + 1), 
                                value=int(np.percentile(self.w, 85) * 100), orientation=QtCore.Qt.Horizontal)
        thr_slider.sliderReleased.connect(Filter)
        thr_slider.valueChanged.connect(Change_Thr)
        thr_slider.setStyleSheet("""
            QSlider{
                background: #E3DEE2;
            }
            QSlider::groove:horizontal {  
                height: 10px;
                margin: 0px;
                border-radius: 5px;
                background: #B0AEB1;
            }
            QSlider::handle:horizontal {
                background: #fff;
                border: 1px solid #E3DEE2;
                width: 17px;
                margin: -5px 0; 
                border-radius: 8px;
            }
            QSlider::sub-page:qlineargradient {
                background: #3B99FC;
                border-radius: 5px;
            }
        """)
        
        label_slider = QW.QLabel()
        
        def Clear():
            if self.Xy['mpl_rects']:
                for rect in self.Xy['mpl_rects']:
                    rect.remove()
            self.Xy = {0: [0, 0, self.hsi.shape[0], self.hsi.shape[1]], 'mpl_rects': [], 'ax_l': 0, 'count': 1} # dct_rect
            combo_box.clear(); combo_box.addItem('Средний спектр изображения', 0); rb_1.setChecked(True)
            Filter()
        
        clear_button = QW.QPushButton('Сбросить ROI')
        clear_button.clicked.connect(Clear)
        Clear()
        
        rb_1.toggled.connect(Draw_Left); rb_2.toggled.connect(Draw_Left)
        mini_h = QW.QHBoxLayout(); mini_h.addWidget(rb_1); mini_h.addWidget(rb_2); mini_h.addWidget(combo_box)
        
        # Комбинироанный виджет
        hbox = QW.QHBoxLayout()
        vbox_1, vbox_2 = QW.QVBoxLayout(), QW.QVBoxLayout()
        vbox_1.addWidget(canvas_l); vbox_1.addWidget(toolbar_l)
        vbox_2.addWidget(clear_button); vbox_2.addWidget(save_button)
        vbox_2.addWidget(self.ui.button_open)
        
        hbox_slider = QW.QHBoxLayout()
        hbox_slider.addWidget(label_slider); hbox_slider.addWidget(thr_slider); hbox_slider.addWidget(QW.QLabel('Маска:'))
        hbox_slider.addWidget(text_w_2); hbox_slider.addWidget(QW.QLabel('-')); hbox_slider.addWidget(text_w_1)
        vbox_2.addLayout(hbox_slider)
        vbox_2.addLayout(mini_h); vbox_2.addWidget(canvas_r)
        
        hbox_down_r = QW.QHBoxLayout(); hbox_down_r.addWidget(QW.QLabel('Границы матрицы:')); 
        hbox_down_r.addWidget(text_1); hbox_down_r.addWidget(text_2)
        hbox_down_r.addWidget(toolbar_r); vbox_2.addLayout(hbox_down_r)
        hbox.addLayout(vbox_1); hbox.addLayout(vbox_2)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QW.QWidget()
        widget.setLayout(hbox)
        self.setCentralWidget(widget)
        self.show()
        
        # Interact
        def onclick(event):
            if str(event.button) == 'MouseButton.RIGHT' and rb_1.isChecked():
                x, y = self.coord(event)
                count = self.Xy['count']
                if count not in self.Xy: self.Xy[count] = []
                self.Xy[count] += [x, y]
                
                if len(self.Xy[count]) > 2:
                    x1, y1, x2, y2 = self.Xy[count]
                    if x1 > x2: x1, x2 = x2, x1
                    if y1 > y2: y1, y2 = y2, y1
                    self.Xy[count] = (x1, y1, x2, y2)
                    patch = ax_l.add_patch(Rectangle(xy=(x1, y1), width=x2-x1, height=y2-y1, 
                                                     fill=False, color='lime', linewidth=2.5))
                    text = ax_l.text((x1+x2)//2, (y1+y2)//2, count, fontsize=15, color='lime')
                    
                    self.Xy['mpl_rects'] += [patch, text]
                    combo_box.addItem(f'ROI {count}', count); combo_box.setCurrentIndex(count)
                    Filter()
                    self.Xy['count'] += 1
                
        fig_l.canvas.mpl_connect('button_press_event', onclick)

if __name__ == "__main__":
    app = QW.QApplication([])
    application = mywindow()
    application.show()
 
    sys.exit(app.exec())