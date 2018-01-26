## From 
## https://github.com/nchah/github-traffic-stats


## OSX
## ------
# pip install --user virtualenv
# mkdir -p ~/my.python.virtualenvs
# cd ~/my.python.virtualenvs
# /Users/admin/Library/Python/2.7/bin/virtualenv gh-traffic-stats  # use any name you please
# ls   # to see the virtual envs
# ##you can switch to another directory on the system if you like
# ##now, the place where you wish to save your traffic stats
source ~/my.python.virtualenvs/gh-traffic-stats/bin/activate # every time you log in to rerun the traffic stats, you will need to rerun this line to 'activate' this virtual env
pip install github_traffic_stats # this will "install" it somewhere inside ~/my.python.virtualenvs/gh-traffic-stats
gts sinarueeger ssimp_software

## HPC
## -------
source ~/my.python.virtualenvs/gh-traffic-stats/bin/activate 
#every time you log in to rerun the traffic stats, you will need to rerun this line to 'activate' this virtual env
pip install github_traffic_stats # this will "install" it somewhere inside ~/my.python.virtualenvs/gh-traffic-stats
gts sinarueeger ssimp_software # This works, you must enter your password

## stores some csv files.

## run B_traffic_github.sh to produce graphs




# ## what aaron wrote on 2018-01-25

#   pip install --user virtualenv
#   # you should consider adding '~/.local/bin' to your $PATH, if you
# haven't already done so
#   # create a directory of your choice to store these virtual envs
#   mkdir -p ~/my.python.virtualenvs
#   cd ~/my.python.virtualenvs
#   # create the virtual environment
#   ~/.local/bin/virtualenv gh-traffic-stats  # use any name you please
#   ls   # to see the virtual envs
#   # you can switch to another directory on the system if you like
# now, the place where you wish to save your traffic stats
#   # source ~/my.python.virtualenvs/gh-traffic-stats/bin/activate #
# every time you log in to rerun the traffic stats, you will need to
# rerun this line to 'activate' this virtual env
#   pip install github_traffic_stats # this will "install" it somewhere
# inside ~/my.python.virtualenvs/gh-traffic-stats
#   gts sinarueeger ssimp_software # This works, you must enter your password

# For me to access your stats, I had to use a slightly different command
# line where I specify the 'organization' via '-o' (i.e. your name as
# you host the repository of interest), then my username (because I know
# my password, not yours!):

#   gts -o sinarueeger aaronmcdaid ssimp_software # i had to run this
# to access your stats

# I guess you should try to run it every week, and then use the merge
# command that is in that traffic stats repo in order "To merge and only
# preserve the unique data points, run: ..."
