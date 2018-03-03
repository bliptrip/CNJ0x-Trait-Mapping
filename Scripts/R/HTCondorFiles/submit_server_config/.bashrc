# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
set -o vi

#Condor command shortcuts
alias cqa='condor_q -analyze'
alias cqh='condor_q -af HoldReason'
alias cr='condor_rm'
alias csi='condor_submit -i'
