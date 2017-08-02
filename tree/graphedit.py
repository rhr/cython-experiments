from matplotlib import patches as Patches
from matplotlib import pyplot as plt

class Draggable(object):
    lock = None
    def __init__(self, ax, label, x=0, y=0, size=20):
        self.ax = ax
        self.label = label
        bbprops = dict(boxstyle='circle', pad=0.5)
        self.text = ax.text(
            x, y, label, ha='center', va='center', size=size, bbox=bbprops)
        self.press = None
        self.background = None

    def connect(self):
        self.cidpress = self.text.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.text.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.text.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.text.axes:
            print('not inaxes')
            return
        if Draggable.lock is not None:
            print('Draggable.lock is not None')
            return
        contains, attrd = self.text.contains(event)
        if not contains:
            print('not contains')
            return
        x0, y0 = self.text.get_position()
        self.press = (x0, y0, event.xdata, event.ydata)

        Draggable.lock = self
        # draw everything but the selected patch and store the pixel buffer
        canvas = self.text.figure.canvas
        axes = self.text.axes
        self.text.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.text)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the patch if the mouse is over us'
        if Draggable.lock is not self:
            return
        if event.inaxes != self.text.axes:
            return
        print(event)
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.text.set_position((x0+dx, y0+dy))

        canvas = self.text.figure.canvas
        axes = self.text.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.text)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        print(event)
        'on release we reset the press data'
        if Draggable.lock is not self:
            return

        self.press = None
        Draggable.lock = None

        # turn off the rect animation property and reset the background
        self.text.set_animated(False)
        self.background = None

        # redraw the full figure
        self.text.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.text.figure.canvas.mpl_disconnect(self.cidpress)
        self.text.figure.canvas.mpl_disconnect(self.cidrelease)
        self.text.figure.canvas.mpl_disconnect(self.cidmotion)

f = plt.figure()
ax = f.add_axes([0,0,1,1])
## e = Patches.Ellipse((0,0), .2, .1)
d = Draggable(ax, label='foo')
d.connect()
