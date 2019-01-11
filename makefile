RSYNC = /usr/local/Cellar/rsync/3.1.3_1/bin/rsync
RSYNCTAGS = --archive --verbose --info=progress2 -au
SRC_FOLDER = $(HOME)/github/sparse_solver/src $(HOME)/github/sparse_solver/julia $(HOME)/github/sparse_solver/main.cpp $(HOME)/github/sparse_solver/tests.cpp

REMOTE = wpq@comps0.cs.toronto.edu

synccode:
	$(RSYNC) $(RSYNCTAGS) $(SRC_FOLDER) $(REMOTE):/u/wpq/github/sparse_solver
