"""
This file contains some low-level useful functions
"""

from __future__ import print_function

from matplotlib.widgets import AxesWidget
import six

def write(*args,**kwargs):
    """
    Wrapper for print function or print to ensure compatibility with python 2
    The arguments are used similarly as the print_function 
    They can also be generalized to python 2 cases
    """
    return print(*args,**kwargs)

#interesting slider for vertical matplotlib
#taken from https://stackoverflow.com/questions/25934279/add-a-vertical-slider-with-matplotlib
class VertSlider(AxesWidget):
    """
    A slider representing a floating point range.

    For the slider to remain responsive you must maintain a
    reference to it.

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*.

        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...).

        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in

        label : str
            Slider label

        valmin : float
            The minimum value of the slider

        valmax : float
            The maximum value of the slider

        valinit : float
            The slider initial position

        label : str
            The slider label

        valfmt : str
            Used to format the slider value, fprint format string

        closedmin : bool
            Indicate whether the slider interval is closed on the bottom

        closedmax : bool
            Indicate whether the slider interval is closed on the top

        slidermin : Slider or None
            Do not allow the current slider to have a value less than
            `slidermin`

        slidermax : Slider or None
            Do not allow the current slider to have a value greater than
            `slidermax`


        dragging : bool
            if the slider can be dragged by the mouse

        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.hline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


def find_files(regexp, archive=None):
    """
    Returns a list of the paths to the files that follow the regular expression
    regexp. They are searched from the current working directory or from an archive 
    given as optional argument.


    :param regexp: A regular expression
    :type regexp: string
    :param archive: an opened tarfile archive (optional)
    :type archive: 
    :returns: a list of all the paths that agree with the regexp
    :rtype: list of strings
    :raises: ValueError if the regexp does not find a single path.


    Example::

        #Find all python files in the current working directory
        find_files('*py') 

        #An exmple outside of the current working directory
        find_files('*/log-*.yaml') 

        #Example using a tarfile
        import tarfile
        my_archive = tarfile.open('archive.tar.gz')
        find_files('*/*/log-*.yaml', archive=my_archive)
    """
    import os

    #Get a list of all paths to files satisfying the regexp
    if archive is not None:
        paths = _find_files_from_archive(regexp, archive) 
    else:
        paths = os.popen('ls '+regexp).read().splitlines()

    #Test that the regexp found files
    if paths == []:
        raise ValueError('The regexp "{}" leads to no file. '\
                         'Consider using another one.'.format(regexp))
    else:
        return paths


def _find_files_from_archive(re, archive):
    """
    This function retrieves the list of Logfiles instances 
    from the file archived satisfying a regular expression.   
    #function to identify an archive out of its regexp, 
    #solves the bug in re for '*' (solved in Python 2.7.6)
    """
    import tarfile    

    #Open the archive
    with tarfile.open(archive, 'r') as arch:
    #Return paths to logfiles satisfying the regexp
        return [f for f in arch.getnames() 
                if all(pattern in f for pattern in re.split('*'))]



if __name__ == '__main__':
    import os

    #Tests of the find_files function
    #
    print("Test finding all python files in this directory")
    print(find_files("*py"))
    print()

    #
    print("Test finding the Utils.py file in this directory")
    print(find_files("Utils.py"))
    print()
 
    #
    print("Test raising a ValueError because the regexp leads to no files")
    try:
        find_files('*html')
    except ValueError as e:
        print('This raised the following ValueError:')
        print(e)
    print()

    #
    print("Test raising an exception because there is no such archive.")
    fname = 'file.tar.gz'
    if fname in os.popen('ls'): os.system('rm '+fname)
    #os.system('rm '+fname)
    try:
        find_files('*py', archive=fname)
    except Exception as e:
        #print(dir(e))
        print('This raised the following Exception:')
        print(e)
    print()

    #
    print("Test without error using an archive")
    os.system('find * -name "*py" | tar -zcvf '+fname+' -T -')
    os.system('ls '+fname)
    find_files('*py', archive=fname)
    os.system('rm '+fname)
