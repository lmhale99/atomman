<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  xmlns="http://www.w3.org/TR/xhtml1/strict">
  <xsl:output method="html" encoding="utf-8" indent="yes" />
  
  <xsl:template match="relaxed-crystal">

    <!-- Shortcut to cell box parameters -->
    <xsl:variable name="ax"><xsl:value-of select="atomic-system/box/avect/value[1]"/></xsl:variable>
    <xsl:variable name="ay"><xsl:value-of select="atomic-system/box/avect/value[2]"/></xsl:variable>
    <xsl:variable name="az"><xsl:value-of select="atomic-system/box/avect/value[3]"/></xsl:variable>
    <xsl:variable name="bx"><xsl:value-of select="atomic-system/box/bvect/value[1]"/></xsl:variable>
    <xsl:variable name="by"><xsl:value-of select="atomic-system/box/bvect/value[2]"/></xsl:variable>
    <xsl:variable name="bz"><xsl:value-of select="atomic-system/box/bvect/value[3]"/></xsl:variable>
    <xsl:variable name="cx"><xsl:value-of select="atomic-system/box/cvect/value[1]"/></xsl:variable>
    <xsl:variable name="cy"><xsl:value-of select="atomic-system/box/cvect/value[2]"/></xsl:variable>
    <xsl:variable name="cz"><xsl:value-of select="atomic-system/box/cvect/value[3]"/></xsl:variable>
    
    <div>
      <!-- Load 3Dmol.js -->
	    <script src='https://3Dmol.csb.pitt.edu/build/3Dmol-min.js'>//</script>

      <style>
        .mol-container {width: 400px; height: 400px; position: relative;}
        .li-spg {list-style-type: circle; list-style-position: inside; margin: 10px;}
        .atomtable {border: 1px solid black; border-collapse: collapse;}
        .atomcol {width: 50px;}
      </style>

      <h1>Relaxed crystal</h1>

      <ul>
        <li><b><xsl:text>UUID4: </xsl:text></b><xsl:value-of select="key"/></li>
        <li><b><xsl:text>Relax method: </xsl:text></b><xsl:value-of select="method"/></li>
        <li><b><xsl:text>Standing: </xsl:text></b><xsl:value-of select="standing"/></li>
        <li><b><xsl:text>Potential: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="potential-LAMMPS/potential/URL">
              <a href="{potential-LAMMPS/potential/URL}"><xsl:value-of select="potential-LAMMPS/potential/id"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="potential-LAMMPS/potential/id"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <li><b><xsl:text>LAMMPS implementation: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="potential-LAMMPS/URL">
              <a href="{potential-LAMMPS/URL}"><xsl:value-of select="potential-LAMMPS/id"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="potential-LAMMPS/id"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <li><b><xsl:text>Family: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="system-info/family-URL">
              <a href="{system-info/family-URL}"><xsl:value-of select="system-info/family"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="system-info/family"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <li><b><xsl:text>Composition: </xsl:text></b><xsl:value-of select="system-info/composition"/></li>
        <li><b><xsl:text>a (Angstrom): </xsl:text></b><xsl:value-of select="system-info/cell/a"/></li>
        <li><b><xsl:text>b (Angstrom): </xsl:text></b><xsl:value-of select="system-info/cell/b"/></li>
        <li><b><xsl:text>c (Angstrom): </xsl:text></b><xsl:value-of select="system-info/cell/c"/></li>
        <li><b><xsl:text>alpha: </xsl:text></b><xsl:value-of select="system-info/cell/alpha"/></li>
        <li><b><xsl:text>beta: </xsl:text></b><xsl:value-of select="system-info/cell/beta"/></li>
        <li><b><xsl:text>gamma: </xsl:text></b><xsl:value-of select="system-info/cell/gamma"/></li>
        <li><b><xsl:text>Epot (eV/atom): </xsl:text></b><xsl:value-of select="potential-energy/value"/></li>
        <li><b><xsl:text>Ecoh (eV/atom): </xsl:text></b><xsl:value-of select="cohesive-energy/value"/></li>
      </ul>

      <p>
        <b><xsl:text>Cell vectors:</xsl:text></b>
        <ul>
          <li>
            <xsl:text>a: [</xsl:text>
            <xsl:value-of select="$ax"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$ay"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$az"/><xsl:text>]</xsl:text>
          </li>
          <li>
            <xsl:text>b: [</xsl:text>
            <xsl:value-of select="$bx"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$by"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$bz"/><xsl:text>]</xsl:text>
          </li>
          <li>
            <xsl:text>c: [</xsl:text>
            <xsl:value-of select="$cx"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$cy"/><xsl:text>, </xsl:text>
            <xsl:value-of select="$cz"/><xsl:text>]</xsl:text>
          </li>
        </ul>
      </p>

      <p>
        <b><xsl:text>Atomic sites:</xsl:text></b>
        <table class="atomtable">
          <tr class="atomtable">
            <th class="atomtable atomcol">id</th>
            <th class="atomtable atomcol">atype</th>
            <th class="atomtable atomcol">x</th>
            <th class="atomtable atomcol">y</th>
            <th class="atomtable atomcol">z</th>
          </tr>
          <xsl:for-each select="atomic-system/atoms/property[name='atype']/data/value">
            <xsl:variable name="row" select="position()"/>
            <xsl:variable name="atype" select="number(current())"/>
            <tr class="atomtable">
              <td class="atomtable atomcol"><xsl:value-of select="$row"/></td>
              <td class="atomtable atomcol"><xsl:value-of select="../../../../../system-info/symbol[$atype]"/></td>
              <!--<td class="atomtable atomcol"><xsl:value-of select="current()"/></td>-->
              <td class="atomtable atomcol"><xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 1], '##0.0000000000000')"/></td>
              <td class="atomtable atomcol"><xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 2], '##0.0000000000000')"/></td>
              <td class="atomtable atomcol"><xsl:value-of select="format-number(../../../property[name='pos']/data/value[($row - 1) * 3 + 3], '##0.0000000000000')"/></td>
            </tr>
          </xsl:for-each>
        </table>

      </p>
    
      <div id="container-01" class="mol-container"><!-- 3D shape goes here! --></div>

      <script>
        
        <!-- Initialize viewer -->
        <xsl:text>$(function() {</xsl:text>
        <xsl:text> let atomcolors = ['', 'blue', 'red', 'green', 'yellow', 'purple', 'black'];</xsl:text>
        <xsl:text> let element = $('#container-01');</xsl:text>
        <xsl:text> let config = { backgroundColor: 'white' };</xsl:text>
        <xsl:text> let viewer = $3Dmol.createViewer( element, config );</xsl:text>

        <!-- Cell line origin to a -->
        <xsl:text> viewer.addLine({start:{x:0.0,y:0.0,z:0.0},end:{x:</xsl:text><xsl:value-of select="$ax"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$ay"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$az"/><xsl:text>}});</xsl:text>

        <!-- Cell line origin to b -->
        <xsl:text> viewer.addLine({start:{x:0.0,y:0.0,z:0.0},end:{x:</xsl:text><xsl:value-of select="$bx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$by"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$bz"/><xsl:text>}});</xsl:text>

        <!-- Cell line origin to c -->
        <xsl:text> viewer.addLine({start:{x:0.0,y:0.0,z:0.0},end:{x:</xsl:text><xsl:value-of select="$cx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$cy"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$cz"/><xsl:text>}});</xsl:text>

        <!-- Cell line a to a+b -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$ax"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$ay"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$az"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $bx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line a to a+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$ax"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$ay"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$az"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $cz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line b to a+b -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$bx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$by"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$bz"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $bx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line b to b+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$bx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$by"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$bz"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($bz + $cz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line c to a+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$cx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$cy"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$cz"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $cz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line c to b+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="$cx"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="$cy"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="$cz"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($bz + $cz)"/><xsl:text>}});</xsl:text>
        
        <!-- Cell line a+b to a+b+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="number($ax + $bx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz)"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz + $cz)"/><xsl:text>}});</xsl:text>
        
        <!-- Cell line a+c to a+b+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="number($ax + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $cz)"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz + $cz)"/><xsl:text>}});</xsl:text>

        <!-- Cell line b+c to a+b+c -->
        <xsl:text> viewer.addLine({start:{x:</xsl:text><xsl:value-of select="number($bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($bz + $cz)"/>
        <xsl:text>},end:{x:</xsl:text><xsl:value-of select="number($ax + $bx + $cx)"/>
        <xsl:text>,y:</xsl:text><xsl:value-of select="number($ay + $by + $cy)"/>
        <xsl:text>,z:</xsl:text><xsl:value-of select="number($az + $bz + $cz)"/><xsl:text>}});</xsl:text>

        <!-- Create atoms -->
        <xsl:for-each select="atomic-system/atoms/property[name='atype']/data/value">
          <xsl:variable name="row" select="position()"/>
          <xsl:text>viewer.addSphere({ center: {x:</xsl:text>
          <xsl:value-of select="number(../../../property[name='pos']/data/value[($row - 1) * 3 + 1])"/>
          <xsl:text>, y:</xsl:text>
          <xsl:value-of select="number(../../../property[name='pos']/data/value[($row - 1) * 3 + 2])"/>
          <xsl:text>, z:</xsl:text>
          <xsl:value-of select="number(../../../property[name='pos']/data/value[($row - 1) * 3 + 3])"/>
          <xsl:text>}, radius: 0.3, color: atomcolors[</xsl:text>
          <xsl:value-of select="."/>
          <xsl:text>] });</xsl:text>
        </xsl:for-each>
        <xsl:text>viewer.zoomTo();viewer.zoom(2);viewer.render();});</xsl:text>
      </script>

    </div>
  </xsl:template>
</xsl:stylesheet>