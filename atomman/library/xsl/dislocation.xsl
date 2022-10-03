<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  xmlns="http://www.w3.org/TR/xhtml1/strict">
  <xsl:output method="html" encoding="utf-8" indent="yes" />
  
  <xsl:template match="dislocation">
    <div>

      <style>
        .aslist {list-style-type: circle; list-style-position: inside; margin: 10px;}
      </style>

      <h1>Dislocation parameter set</h1>
      <ul>
        <li><b><xsl:text>ID: </xsl:text></b><xsl:value-of select="id"/></li>
        <li><b><xsl:text>UUID4: </xsl:text></b><xsl:value-of select="key"/></li>
        <li><b><xsl:text>Family: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="system-family-URL">
              <a href="{system-family-URL}"><xsl:value-of select="system-family"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="system-family"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <li><b><xsl:text>Burgers vector: </xsl:text></b><xsl:value-of select="Burgers-vector"/></li>
        <li>
          <b><xsl:text>Slip plane: </xsl:text></b><xsl:text>(</xsl:text>
          <xsl:for-each select="slip-plane">
            <xsl:value-of select="."/>
            <xsl:if test="position() &lt; last()">
              <xsl:text>, </xsl:text>
            </xsl:if>
          </xsl:for-each>
          <xsl:text>)</xsl:text>
        </li>
        <li>
          <b><xsl:text>Line direction: </xsl:text></b><xsl:text>[</xsl:text>
          <xsl:for-each select="line-direction">
            <xsl:value-of select="."/>
            <xsl:if test="position() &lt; last()">
              <xsl:text>, </xsl:text>
            </xsl:if>
          </xsl:for-each>
          <xsl:text>]</xsl:text>
        </li>
        <li><b><xsl:text>Character: </xsl:text></b><xsl:value-of select="character"/></li>
        <li>
          <b><xsl:text>atomman Dislocation defect parameters:</xsl:text></b>
          <ul class="aslist">
            <li><xsl:text>slip_hkl: "</xsl:text><xsl:value-of select="calculation-parameter/slip_hkl"/><xsl:text>"</xsl:text></li>
            <li><xsl:text>ξ_uvw: "</xsl:text><xsl:value-of select="calculation-parameter/ξ_uvw"/><xsl:text>"</xsl:text></li>
            <li><xsl:text>burgers: "</xsl:text><xsl:value-of select="calculation-parameter/burgers"/><xsl:text>"</xsl:text></li>
            <li><xsl:text>m: "</xsl:text><xsl:value-of select="calculation-parameter/m"/><xsl:text>"</xsl:text></li>
            <li><xsl:text>n: "</xsl:text><xsl:value-of select="calculation-parameter/n"/><xsl:text>"</xsl:text></li>
            <xsl:if test="calculation-parameter/shift">
              <li><xsl:text>shift: "</xsl:text><xsl:value-of select="calculation-parameter/shift"/><xsl:text>"</xsl:text></li>
            </xsl:if>
            <xsl:if test="calculation-parameter/shiftscale">
              <li><xsl:text>shiftscale: "</xsl:text><xsl:value-of select="calculation-parameter/shiftscale"/><xsl:text>"</xsl:text></li>
            </xsl:if>
            <xsl:if test="calculation-parameter/shiftindex">
              <li><xsl:text>shiftindex: "</xsl:text><xsl:value-of select="calculation-parameter/shiftindex"/><xsl:text>"</xsl:text></li>
            </xsl:if>
          </ul>
        </li>
      </ul>
      
    </div>

  </xsl:template>
</xsl:stylesheet>