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
        self.bbpatch = self.text.get_bbox_patch()
        self.press = None
        self.background = None

    def connect(self):
        self.cidpress = self.patch.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.patch.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.patch.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.patch.axes:
            return
        if Draggable.lock is not None:
            return
        contains, attrd = self.patch.contains(event)
        if not contains:
            return
        x0, y0 = self.patch.center
        self.press = x0, y0, event.xdata, event.ydata

        Draggable.lock = self
        # draw everything but the selected patch and store the pixel buffer
        canvas = self.patch.figure.canvas
        axes = self.patch.axes
        self.patch.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.patch)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the patch if the mouse is over us'
        if Draggable.lock is not self:
            return
        if event.inaxes != self.patch.axes:
            return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.patch.center = (x0+dx, y0+dy)

        canvas = self.patch.figure.canvas
        axes = self.patch.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.patch)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if Draggable.lock is not self:
            return

        self.press = None
        Draggable.lock = None

        # turn off the rect animation property and reset the background
        self.patch.set_animated(False)
        self.background = None

        # redraw the full figure
        self.patch.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.patch.figure.canvas.mpl_disconnect(self.cidpress)
        self.patch.figure.canvas.mpl_disconnect(self.cidrelease)
        self.patch.figure.canvas.mpl_disconnect(self.cidmotion)

f = plt.figure()
ax = f.add_axes([0,0,1,1])
e = Patches.Ellipse((0,0), .2, .1)
d = Draggable(e)
