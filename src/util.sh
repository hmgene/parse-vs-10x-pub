       
repl()  
{   
    usage="$FUNCNAME <file> <key_value.txt> [replaced-only=0]
    multipe matches are delimitaated by |
    ";
    if [ $# -lt 2 ]; then echo "$usage"; return; fi;
    local tmp=$(mktemp -d );
    cat $2 > $tmp/a;
    cat $1 | perl -e 'use strict; my $only='${3:-0}';
                open(my $fh,"<","'$tmp/a'") or die $!;
        my %h=();map{chomp;my($k,$v)=split/\t/,$_; $h{$k}{$v}++;} <$fh>;
                close($fh);
                while(<STDIN>){chomp;  my @d=split/\t/,$_;
            my $hit=0;
            my @r=map {
                if( defined $h{$_} ){
                    $hit++;
                    $_=join("|",keys %{$h{$_}});
                } 
                $_ 
            } @d;
            if( $hit > 0 || $only == 0 ){
                print join("\t",@r),"\n";
            }

                }
        '
}

