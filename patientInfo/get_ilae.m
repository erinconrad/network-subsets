function outcome = get_ilae(name)

all_names = [{'HUP064'  }
    {'HUP065'  }
    {'HUP068'  }
    {'HUP070'  }
    {'HUP073'  }
    {'HUP074'  }
    {'HUP075'  }
    {'HUP078'  }
    {'HUP080'  }
    {'HUP082'  }
    {'HUP083'  }
    {'HUP086'  }
    {'HUP087'  }
    {'HUP088'  }
    {'HUP094'  }
    {'HUP105'  }
    {'HUP106'  }
    {'HUP107'  }
    {'HUP111A' }
    {'HUP111B' }
    {'HUP116'  }
    {'Study012'}
    {'Study016'}
    {'Study017'}
    {'Study019'}
    {'Study020'}
    {'Study022'}
    {'Study028'}
    {'Study029'}];


all_outcomes = [1
1
1
2
1
1
5
4
2
1
2
1
3
1
2
1
2
1
1
1
1
1
4
4
5
4
5
4
5];

found_it = 0;
for i = 1:length(all_names)
    if strcmp(name,all_names{i}) == 1
        found_it = 1;
        index = i;
        break
    end
end

if found_it == 0
    fprintf('Did not find patient %s.\n',name);
    outcome = nan;
    return 
end

outcome = all_outcomes(index);

end