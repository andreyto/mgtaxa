#!/bin/sh

if [ $1 = -h ]; then 

	echo '

 		This to work needs in your system to have installed mercurial, git and additionally:

 		-- dulwitch 0.6.0 (released 05/22/10) or later version, which contains Checksum bugfix. 
    	       	   Get it from here https://launchpad.net/dulwich

 		-- Mercurial to Git plugin, http://hg-git.github.com/
    		   Can be installed using Python easy install toolkit and running:  easy_install hg-git


 		TO TEST THE MERCURIAL TO GITHUB PROCESS, WE FOLLOW THE STEPS BELOW

 		-- Fork django piston http://bitbucket.org/jespern/django-piston/wiki/Home  ** as hg-django-piston **
    		   This will be used instead of Galaxy for the tests since it is small

 		-- Create an empty Github repository called  ** django-piston **

 		-- Finally execute this script again with your github username

		Options:

		-h : prints this help menu
		-u github_user : your github username to start the repository transfer process

		'

fi


if [ $1 = -u ]; then

        echo
	echo
	echo '  *** Here I am Ntino and I clone the Galaxy repo  *** '
        echo
	hg clone https://$2@bitbucket.org/$2/hg-django-piston
        cd hg-django-piston
	hg bookmark -r default master



        echo
        echo
	echo
	echo '  ***  Push galaxy from mercurial to github. The mercurial to github plugin (hg-git)'
	echo 'takes care of the transfer  *** '
        echo
	hg push git+ssh://git@github.com:$2/django-piston.git

	
        echo
        echo
	echo
	echo '  ***  Here the Galaxy developers add some files   *** '
        echo
	cd ..
	mkdir devs-hg-django-piston/
	cd devs-hg-django-piston/
	hg clone https://$2@bitbucket.org/$2/hg-django-piston
	cd hg-django-piston/
	touch test-mercurial-to-git-1
	hg add test-mercurial-to-git-1
	hg commit -m 'test-mercurial-git-1'
	hg push


        echo
        echo
        echo
        echo
        echo '  ***  Here I am Ntino again (going back to Ntino Galaxy clone), add some files, and send them to my version of Galaxy on Github  *** '
        echo
        cd ../../hg-django-piston
	touch test-git-to-mercurial-1
	hg add test-git-to-mercurial-1
	hg commit -m 'test-git-to-mercurial-1'
	hg push git+ssh://git@github.com:$2/django-piston.git


        echo
	echo
	echo '  ***  Here we merge the changes that the Galaxy developers added with the changes Ntino did on Github  *** '
 	hg pull
	hg merge
	hg commit -m 'merging changes that Galaxy developers did with my code on github'
	hg push git+ssh://git@github.com:$2/django-piston.git		


	echo
	echo ' All works !!! '
	echo
	echo

fi

exit
