#!/usr/bin/env bash

patch_branch=0

# command-line options
while [[ -n $1 ]]; do
	case "$1" in
		--branch|branch)
			patch_branch=1
		;;
		*)
			break
		;;
	esac
	shift
done

# positional args
if [[ -n $1 ]]; then
	alltut=( "${@}" )
else
	alltut=( ./tutorial*.md )
fi

echo "Will run on:"
printf "<%s>\n" "${alltut[@]}"

# branch (for usual PR and pushes)
branch="$(git branch --show-current)"
# fall back on tags if branch was left empty (checked out in grafted/detached head mode, e.g. when testing a tag)
branch="${branch:-$(git describe --tags)}"
# fall back to hash
branch="${branch:-$(git rev-parse HEAD)}"

default="master"
echo "current branch: ${branch}"


if (( patch_branch )); then
	if [[ -z "${branch}" ]]; then
		echo -e "\e[31;1mCannot determine current branch or tag!\e[0m The installer in tutorials will still clone the \`${default}\` branch."
	elif [[ "${branch}" != "${default}" ]]; then
		echo "patching tutorials..."
		sed -ri "s@(quick_install.sh) +(-[^b])@\1 -b ${branch} \2@g;s@(V-pipe)/${default}(/utils)@\1/${branch}\2@g" "${alltut[@]}"
		## example command to keep the output with errors
		# jupytext --to ipynb --pipe-fmt ipynb --pipe 'jupyter nbconvert --to ipynb --execute --allow-errors --stdin --stdout' ./*.md
	else
		echo -e "\e[32;1mNo need to patch branch.\e[0m The default branch is already \`${branch}\`"
	fi
elif [[ "${branch}" != "${default}" ]]; then
	echo -e "\e[33;1mWarning: installer in tutorials will still clone the \`${default}\` branch!\e[0m Use \`$0 branch\` to patch them on the fly and instll the current branch \`${branch}\`"
fi

# create Jupyter Notebooks for all Markdown files
exec jupytext --set-formats ipynb,md --execute "${alltut[@]}"
