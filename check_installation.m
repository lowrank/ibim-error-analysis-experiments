function ret = check_installation()
    [~, cmdout] = system('if command -v wolframscript > /dev/null 2>&1; then echo "1"; else echo "0"; fi');
    ret = str2double(cmdout) == 1;
end
