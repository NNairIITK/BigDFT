# jhbuild - a tool to ease building collections of source packages
# Copyright (C) 2001-2006  James Henstridge
#
#   base.py: the most common jhbuild commands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from optparse import make_option

import jhbuild.moduleset
import jhbuild.frontends
from jhbuild.errors import UsageError, FatalError
from jhbuild.commands import Command, register_command


class cmd_dist(Command):
    doc = _('Create tar packages for all modules')

    name = 'dist'
    usage_args = '[ options ... ] [ modules ... ]'

    def __init__(self):
        Command.__init__(self, [
            make_option('-s', '--skip', metavar='MODULES',
                        action='append', dest='skip', default=[],
                        help=_('don\'t package the given modules')),
            make_option('-t', '--start-at', metavar='MODULE',
                        action='store', dest='startat', default=None,
                        help=_('start building at the given module')),
            make_option('-c', '--compress-kind', metavar='COMP',
                        action='store', dest='compress', default="gz",
                        help=_('compression type (default is gzip)')),
            make_option('-f', '--full',
                        action='store_true', dest='full', default=False,
                        help=_('package also non BigDFT specific libraries')),
            ])

    def run(self, config, options, args, help=None):
        for item in options.skip:
            config.skip += item.split(',')

        module_set = jhbuild.moduleset.load(config)
        module_list = module_set.get_module_list(args or config.modules, config.skip)
        # remove modules up to startat
        if options.startat:
            while module_list and module_list[0].name != options.startat:
                del module_list[0]
            if not module_list:
                raise FatalError(_('%s not in module list') % options.startat)

        config.compress = options.compress
        config.full = options.full

        build = jhbuild.frontends.get_buildscript(config, module_list, module_set=module_set)
        return build.build(phases=['dist'])

register_command(cmd_dist)
