## Note su Worhp
scaricare dal sito da https://worhp.de/latest/ la versione per Ubuntu 20.04 (`worhp_1.14-0_ubuntu2004.deb`)

installare con il comando

    sudo dpkg -i worhp_1.14-0_ubuntu2004.deb

In case Ubuntu complains about missing dependance, run

    sudo apt --fix-broken install
 
Se stai usando WSL2, ad ogni reboot cambierà il MAC address e quindi non funzionerà

Puoi forzare un MAC address `XX:XX:XX:XX:XX:XX` arbitrario in wsl2 con il comando 

    sudo ip link set dev bond0 address XX:XX:XX:XX:XX:XX
