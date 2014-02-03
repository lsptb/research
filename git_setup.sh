#!/bin/bash

# check git repo path argument
dot_dir=$1
if [ -z "${dot_dir}" ]; then
    echo "No argument supplied. Try from the directory: sh git_setup.sh \`pwd\`"
elif [ ! -d "${dot_dir}" ]; then
    echo "Argument is not a directory. Try from the directory: sh git_setup.sh \`pwd\`"
else
    # if everything looks ok...
    # softlink gitconfig into home directory, will not overwrite files
    file="gitconfig"
    backup="${HOME}/.gitconfig_old"
    # check if the file is already symlinked
    if  [ -h "${HOME}/.${file}" ]; then
	echo "${HOME}/.${file} is already symlinked:"
	ls -al ${HOME}/.${file}
    # check if the file exists and move before linking
    elif [ -f "${HOME}/.${file}" ]; then
	mv ${HOME}/.${file} ${backup}
	echo "${HOME}/.${file} already exists. It has been moved to ${backup}"
	ln -s $dot_dir/$file ${HOME}/.$file 2> /dev/null
	echo "File link created:"
	ls -al ${HOME}/.${file}
    # otherwise, symlink it
    else
	ln -s $dot_dir/$file ${HOME}/.$file 2> /dev/null
	echo "File link created:"
	ls -al ${HOME}/.${file}
    fi
    # set gitignore global to the versioned file in the dotfiles directory
    git config --global core.excludesfile ${dot_dir}/gitignore_global

    # get user's name and email
    echo "Enter a name to appear on your git commits: "
    read name
    echo "Enter an email address for your git commits: "
    read email

    # add name and email to git config
    git config --global user.name "${name}"
    git config --global user.email "${email}"
    git config --global core.editor vim
fi
