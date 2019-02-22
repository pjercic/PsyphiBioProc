# PsyphiBioProc
Psychophysiology Biofeedback Processor

==============================
 What is psyphibio-proc?
==============================

psyphibio-proc is a complete psychophisiology processor.

At this time, support for C++, linux, windows.

You have to use Make in linux and VS2010 in windows.

=================
 How To Get gogo
=================

Download the latest version from
http://bitbucket.org/pedersen/gogo/downloads (you may choose a
released version, or use the current latest in development version,
named "tip".

==============================
 Configuring Your Environment
==============================

gogo uses four environment variables:

 * EDITOR : You should likely have this already set. If not, make sure
   to set it now, and point it at your favorite editor. I recommend
   EMACS
 * PRJDIRS : This is a space separated list of directories in which
   projects can be stored. Default: ${HOME}/projects
 * NEWPRJ : This is the directory in which gogo will create new
   project environments. Default: ${HOME}/projects
 * GOGO_PLUGIN_DIR : This is the directory where gogo plugins are
   located. Plugin files are simply bash shell scripts that have a
   name ending with ".gogo". They will be sourced into the running
   shell. Default: <location_of_.bash.gogo>/gogoplugins

If you are okay with these variables, you don't have to do
anything. Otherwise, edit ${HOME}/.bashrc and set the values to your
liking. Make sure to restart bash before continuing.

=================
 Installing gogo
=================

Extract the download into a directory. I recommend doing this in
${NEWPRJ}, so that you will have a folder named ${NEWPRJ}/gogo Keep
all the files around, it will be useful (and the total size of the
files is quite small). In ${HOME}/.bashrc, add the following line:

  source ${NEWPRJ}/gogo/.bash.gogo

That is all that you*must* do to make gogo available to you. Of
course, it won't be terribly useful to you until you get other tools
installed. The recommended tool chain at this time:

 * Mercurial : http://mercurial.selenic.com/
 * Python : http://www.python.org/
 * virtualenv : http://virtualenv.openplans.org/
 * virtualenvwrapper : http://www.doughellmann.com/projects/virtualenvwrapper/
 * exuberant ctags : http://ctags.sourceforge.net/
 * emacs : http://www.gnu.org/software/emacs/
 *  running in daemon mode : http://www.gnu.org/software/emacs/manual/html_node/emacs/Emacs-Server.html#Emacs-Server

You can run "gogo -k" to check your environment, and have it report on
any missing items that could affect you.

===============
 EMACS Support
===============

If you are using EMACS, and are running it in daemon mode, you have a
couple of extra options available to you. I recommend making sure to
load desktops.el and virtualenv.el in your .emacs file. These two
files are key in the usage of gogo with EMACS.

When you do use these, you will find that your inferior Python
interpreter is automatically set to use the python executable found in
your virtualenv for the project, and this will be updated as you
switch from project to project that uses a Python virtualenv.

============
 gogo Usage
============

You have three commands available to you:

* gogo
* gogo_build_tags
* gcd

gogo
====

Once you have completed all the installation steps, run "gogo -h" to
see what the usage and options are. To provide an example of how to
use it, I'll show you how to clone gogo from the original bitbucket
repository, and update it easily.

  gogo -n -c http://bitbucket.org/pedersen/gogo gogo

This command will create a new project (-n), clone the bitbucket
repository (-c url), and put the results into ${NEWPRJ}/gogo.

  gogo -u gogo

This will download updates from the gogo bitbucket repository.

  gogo -d gogo

This will ask for confirmation, and (if you say "y"), destroy the gogo
checkout.

A longer example is available thanks to the TurboGears project
( http://www.turbogears.org/ ). Building the documentation manually
takes a bit of effort. You have to create a new virtualenv, install
the required tools, install the required packages for the
documentation, and then you can finally edit the files themselves.

Using gogo, the process couldn't be simpler:

  gogo -n -l python -t tg21docs tg21docs

This will create a new project (-n), create a python virtualenv (-l
python), create a project using the tg21docs project template (-t
tg21docs), clone the documentation repository, and put that clone into
${NEWPRJ}/tg21docs. Finally, it will activate the virtualenv, and cd
into the top directory of the project. If you are using EMACS in
daemon mode, EMACS will even be updated to be looking at the files
(and version control) for the project.

When you come back in a week, and need to pull down new updates, it's
still very simple:

  gogo -u tg21docs

gogo will activate the virtualenv, run "hg pull -u" for you, and cd
into ${NEWPRJ}/tg21docs.

Finally, when it comes time to delete:

  gogo -d -f tg21docs

This will destroy the tg21docs virtualenv *and*
${NEWPRJ}/tg21docs. Your system will be as if the project was never
installed.

gogo_build_tags
===============

You must activate a project before using this command. Once you have
done that, though, it does not matter where in your file system you
are. Running this command will rebuild the tags file for your project.

gcd
===

You must activate a project before using this command. Once you have
done that, though, it does not matter where in your file system you
are. Using gcd by itself will cd to the top level of that
project. Passing in a path relative to that top path will cd to that
directory. So, if you have <project>/my/deep/dir/here, after
activating <project>, you may then use "gogo my/deep" and your prompt
will now be at <project>/my/deep regardless of where you were in the
file system.

=====================
 Additional Features
=====================

 * Individual plugins can be enabled and disabled. This is
   accomplished with the command gogo_enable and gogo_disable. Note
   that plugins immediate in ${GOGO_PLUGIN_DIR} are considered core,
   and cannot be disabled. Any plugins placed in subdirectories under
   ${GOGO_PLUGIN_DIR}, though, can be disabled and enabled. You will
   simply need to restart your shell for it to take effect.
 * You may view the list of plugins, and see which ones are enabled
   and disabled, with the command gogo_list_plugins.
 * gogo, gcd, gogo_enabled, and gogo_disable all support tab
   completion. When you hit tab, gogo will try to complete with
   project names.
 * gogo will automatically rebuild a tags file for you. Whenever you
   switch to a project, gogo automatically runs ctags before
   returning. Furthermore, you can always re-run the tag building
   after switching, just by using the command: gogo_build_tags
 * If you have a <project>/.gogo file, then that file will be sourced
   when the project is activated. It is also recommended to have a
   project_deactivate in this file which will undo the environmental
   changes caused by it. For instance, deleting any user defined
   functions you create. This function will be called when the project
   is deactivated. Note that you do not need to delete
   project_deactivate itself; gogo will do that for you.

=======================
 Inspiration & Credits
=======================

I have to admit that inspiration and code came from external sources.

 * desktops.el : This mostly came from
   http://stackoverflow.com/questions/847962/what-alternate-session-managers-are-available-for-emacs/849180#849180
   (further information at
   http://scottfrazersblog.blogspot.com/2009/12/emacs-named-desktop-sessions.html
   ), and I added some very minor enhancements to remember the name of
   the currently loaded desktop.
 * virtualenv.el : This mostly came from
   http://jesselegg.com/archives/2010/03/14/emacs-python-programmers-2-virtualenv-ipython-daemon-mode/
   and was minorly enhanced by me to allow the virtualenv-activate
   functions to take optional virtualenv names

Believe it or not, that's where my inspiration came from.

===========
 ToDo List
===========

 * Add support for emacs (not daemon mode), vim, and vi (yes, separately)
 * Add support for other VCS systems (git, svn)
 * Find way to learn where the desktop-sessions directory is from .emacs


# Local Variables:  #
# mode: rst         #
# End:              #
